# Reimplementation of BuildHON (Xu, Wickramarathne and Chawla 2016; Saebi et al. 2020) written for Nestimate's
# equivalence tests. Moved verbatim from local_testing_and_equivalence/test-hon-pyhon-equiv.R. Standard library only.
from collections import defaultdict, Counter
import math

def run_buildhon(trajectories, max_order, min_support):
    # --- Stage 1: Extract Rules ---
    Count = {}
    Distribution = defaultdict(dict)
    SourceToExtSource = {}
    Rules = defaultdict(dict)

    # BuildObservations
    for traj in trajectories:
        n = len(traj)
        for order in range(2, max_order + 2):
            if order > n:
                break
            for start in range(n - order + 1):
                subseq = tuple(traj[start:start + order])
                Target = subseq[-1]
                Source = subseq[:-1]
                if Source not in Count:
                    Count[Source] = Counter()
                Count[Source][Target] += 1

    # BuildDistributions
    for Source in Count:
        for Target in list(Count[Source].keys()):
            if Count[Source][Target] < min_support:
                Count[Source][Target] = 0
        for Target in Count[Source]:
            if Count[Source][Target] > 0:
                Distribution[Source][Target] = Count[Source][Target] / sum(Count[Source].values())

    # BuildSourceToExtSource
    for source in Distribution:
        if len(source) > 1:
            NewOrder = len(source)
            for starting in range(1, len(source)):
                curr = source[starting:]
                if curr not in SourceToExtSource:
                    SourceToExtSource[curr] = {}
                if NewOrder not in SourceToExtSource[curr]:
                    SourceToExtSource[curr][NewOrder] = set()
                SourceToExtSource[curr][NewOrder].add(source)

    def ExtendSource(Curr, NewOrder):
        if Curr in SourceToExtSource:
            if NewOrder in SourceToExtSource[Curr]:
                return SourceToExtSource[Curr][NewOrder]
        return []

    def KLD(a, b):
        divergence = 0
        for target in a:
            pa = a.get(target, 0)
            pb = b.get(target, 0)
            if pa > 0 and pb > 0:
                divergence += pa * math.log(pa / pb, 2)
            elif pa > 0 and pb == 0:
                divergence = float('inf')
        return divergence

    def KLDThreshold(NewOrder, ExtSource):
        return NewOrder / math.log(1 + sum(Count[ExtSource].values()), 2)

    def AddToRules(Source):
        if len(Source) > 0:
            Rules[Source] = dict(Distribution[Source])
            PrevSource = Source[:-1]
            AddToRules(PrevSource)

    def ExtendRule(Valid, Curr, order, MaxOrder):
        if order >= MaxOrder:
            AddToRules(Valid)
        else:
            Distr = Distribution[Valid]
            NewOrder = order + 1
            Extended = ExtendSource(Curr, NewOrder)
            if len(Extended) == 0:
                AddToRules(Valid)
            else:
                for ExtSource in Extended:
                    ExtDistr = Distribution[ExtSource]
                    if KLD(ExtDistr, Distr) > KLDThreshold(NewOrder, ExtSource):
                        ExtendRule(ExtSource, ExtSource, NewOrder, MaxOrder)
                    else:
                        ExtendRule(Valid, ExtSource, NewOrder, MaxOrder)

    # GenerateAllRules
    for Source in list(Distribution.keys()):
        if len(Source) == 1:
            AddToRules(Source)
            ExtendRule(Source, Source, 1, max_order)

    # --- Stage 2: Build Network ---
    Graph = defaultdict(dict)
    SortedSource = sorted(Rules, key=lambda x: len(x))
    for source in SortedSource:
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                PrevSource = source[:-1]
                PrevTarget = (source[-1],)
                if PrevSource not in Graph or source not in Graph[PrevSource]:
                    if PrevTarget in Graph.get(PrevSource, {}):
                        Graph[PrevSource][source] = Graph[PrevSource][PrevTarget]
                        del Graph[PrevSource][PrevTarget]

    # RewireTails
    ToAdd = []
    ToRemove = []
    for source in list(Graph.keys()):
        for target in list(Graph[source].keys()):
            if len(target) == 1:
                NewTarget = source + target
                while len(NewTarget) > 1:
                    if NewTarget in Graph:
                        ToAdd.append((source, NewTarget, Graph[source][target]))
                        ToRemove.append((source, target))
                        break
                    else:
                        NewTarget = NewTarget[1:]
    for (source, target, weight) in ToAdd:
        Graph[source][target] = weight
    for (source, target) in ToRemove:
        if target in Graph.get(source, {}):
            del Graph[source][target]

    # Convert to edge list
    def SequenceToNode(seq):
        curr = seq[-1]
        node = curr + '|'
        seq = seq[:-1]
        while len(seq) > 0:
            curr = seq[-1]
            node = node + curr + '.'
            seq = seq[:-1]
        if node[-1] == '.':
            return node[:-1]
        else:
            return node

    edges = []
    for source in Graph:
        for target in Graph[source]:
            edges.append({
                'from': SequenceToNode(source),
                'to': SequenceToNode(target),
                'weight': Graph[source][target]
            })

    return {'edges': edges, 'n_edges': len(edges)}
