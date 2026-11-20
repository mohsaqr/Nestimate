# Reimplementation of BuildHON+ (Xu, Wickramarathne and Chawla 2016; Saebi et al. 2020) written for Nestimate's
# equivalence tests. Moved verbatim from local_testing_and_equivalence/test-hon-pyhon-equiv.R. Standard library only.
import math
from collections import defaultdict

def run_hon_plus(trajectories, max_order, min_freq):
    Count = defaultdict(lambda: defaultdict(int))
    Distribution = defaultdict(dict)
    Rules = defaultdict(dict)
    SourceToExtSource = defaultdict(set)
    StartingPoints = defaultdict(set)

    # Build order 1
    for ti in range(len(trajectories)):
        traj = trajectories[ti]
        for index in range(len(traj) - 1):
            Source = tuple(traj[index:index+1])
            Target = traj[index+1]
            Count[Source][Target] += 1
            StartingPoints[Source].add((ti, index))

    for Source in list(Count.keys()):
        if len(Source) == 1:
            for Target in list(Count[Source].keys()):
                if Count[Source][Target] < min_freq:
                    Count[Source][Target] = 0
            total = sum(Count[Source].values())
            if total > 0:
                for Target in Count[Source]:
                    if Count[Source][Target] > 0:
                        Distribution[Source][Target] = 1.0 * Count[Source][Target] / total

    def KLD(a, b):
        d = 0
        for t in a:
            pa = a.get(t, 0)
            pb = b.get(t, 0)
            if pa > 0 and pb > 0:
                d += pa * math.log(pa / pb, 2)
            elif pa > 0:
                d = float('inf')
        return d

    def KLDThreshold(no, es):
        total = sum(Count[es].values())
        if total == 0:
            return float('inf')
        return no / math.log(1 + total, 2)

    def MaxDivergence(Distr):
        mk = sorted(Distr, key=Distr.__getitem__)
        return {mk[0]: 1}

    def ExtendObservation(Source):
        if len(Source) > 1:
            if Source[1:] not in Count or len(Count[Source]) == 0:
                ExtendObservation(Source[1:])
        order = len(Source)
        C = defaultdict(lambda: defaultdict(int))
        for ti, index in StartingPoints[Source]:
            traj = trajectories[ti]
            if index - 1 >= 0 and index + order < len(traj):
                ExtSource = tuple(traj[index-1:index+order])
                Target = traj[index+order]
                C[ExtSource][Target] += 1
                StartingPoints[ExtSource].add((ti, index - 1))
        if len(C) == 0:
            return
        for s in C:
            for t in C[s]:
                if C[s][t] < min_freq:
                    C[s][t] = 0
                Count[s][t] += C[s][t]
            total = sum(C[s].values())
            if total > 0:
                for t in C[s]:
                    if C[s][t] > 0:
                        Distribution[s][t] = 1.0 * C[s][t] / total
                        SourceToExtSource[s[1:]].add(s)

    def ExtendSourceFast(Curr):
        if Curr in SourceToExtSource:
            return SourceToExtSource[Curr]
        else:
            ExtendObservation(Curr)
            if Curr in SourceToExtSource:
                return SourceToExtSource[Curr]
            else:
                return []

    def AddToRules(Source):
        for order in range(1, len(Source)+1):
            s = Source[0:order]
            if s not in Distribution or len(Distribution[s]) == 0:
                ExtendSourceFast(s[1:])
            Rules[s] = Distribution[s]

    def ExtendRule(Valid, Curr, order):
        if order >= max_order:
            AddToRules(Valid)
        else:
            Distr = Distribution[Valid]
            CurrDistr = Distribution.get(Curr, {})
            if len(CurrDistr) == 0:
                AddToRules(Valid)
                return
            if KLD(MaxDivergence(CurrDistr), Distr) < KLDThreshold(order+1, Curr):
                AddToRules(Valid)
            else:
                Extended = ExtendSourceFast(Curr)
                if len(Extended) == 0:
                    AddToRules(Valid)
                else:
                    for ext in Extended:
                        ed = Distribution.get(ext, {})
                        if len(ed) > 0 and KLD(ed, Distr) > KLDThreshold(order+1, ext):
                            ExtendRule(ext, ext, order+1)
                        else:
                            ExtendRule(Valid, ext, order+1)

    for Source in tuple(Distribution.keys()):
        AddToRules(Source)
        ExtendRule(Source, Source, 1)

    # Build network (same as BuildNetwork.py)
    Graph = defaultdict(dict)

    def SequenceToNode(seq):
        l = seq[-1]
        if len(seq) == 1:
            return l + '|'
        return l + '|' + '.'.join(reversed(seq[:-1]))

    def Rewire(source, target):
        ps = source[:-1]
        pt = (source[-1],)
        if ps in Graph and source not in Graph[ps]:
            if pt in Graph.get(ps, {}):
                Graph[ps][source] = Graph[ps][pt]
                del Graph[ps][pt]

    SortedSource = sorted(Rules, key=lambda x: len(x))
    for source in SortedSource:
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                Rewire(source, (target,))

    # Tail rewiring
    ta = []
    tr = []
    for source in Graph:
        for target in list(Graph[source].keys()):
            if len(target) == 1:
                nt = source + target
                while len(nt) > 1:
                    if nt in Graph:
                        ta.append((source, nt, Graph[source][target]))
                        tr.append((source, target))
                        break
                    nt = nt[1:]
    for s, t, w in ta:
        Graph[s][t] = w
    for s, t in tr:
        if t in Graph.get(s, {}):
            del Graph[s][t]

    return [{'from': SequenceToNode(s), 'to': SequenceToNode(t), 'weight': Graph[s][t]}
            for s in Graph for t in Graph[s]]
