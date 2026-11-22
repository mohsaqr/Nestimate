# Runs the authors' BuildHON (reference/buildhon/BuildRules.py, Xu, Wickramarathne and Chawla 2016) or BuildHON+
# (reference/buildhon/BuildRulesFastParameterFree.py, Saebi et al. 2020) and their BuildNetwork.py, as vendored
# (unchanged, commit 2c89ad4 of github.com/xyjprc/hon). Written for Nestimate's equivalence tests.
# The authors' modules keep their state in module globals (BuildRules.Initialize() and BuildNetwork.Initialize()
# do not reset them), so every run reloads them; BuildHON+ may write divergences.csv to the working directory, so
# every run works in a fresh temporary directory. Higher-order nodes are written "A -> B" (oldest state first).
import importlib
import os
import sys
import tempfile

def _node(sequence):
    return ' -> '.join(sequence)

def run_buildhon_reference(folder, trajectories, max_order, min_support, plus):
    folder = os.path.abspath(folder)
    if folder not in sys.path:
        sys.path.insert(0, folder)
    import BuildRules, BuildRulesFastParameterFree, BuildNetwork
    rules_module = importlib.reload(BuildRulesFastParameterFree if plus else BuildRules)
    network_module = importlib.reload(BuildNetwork)
    rules_module.Verbose = False
    network_module.Verbose = False
    records = [[str(i), [str(s) for s in t]] for i, t in enumerate(trajectories)]
    here = os.getcwd()
    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        try:
            rules = rules_module.ExtractRules(records, int(max_order), int(min_support))
            # BuildNetwork() changes the rules it reads ("destructive to Rules"): read them first.
            rule_rows = [{'from': _node(source), 'to': str(target), 'probability': float(p)}
                         for source in rules for target, p in rules[source].items()]
            network = network_module.BuildNetwork(rules)
        finally:
            os.chdir(here)
    edge_rows = [{'from': _node(source), 'to': _node(target), 'weight': float(w)}
                 for source in network for target, w in network[source].items()]
    return {'rules': rule_rows, 'edges': edge_rows}
