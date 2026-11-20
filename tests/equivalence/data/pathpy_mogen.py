# pathpy 2.2.0 (Scholtes 2017) read for Nestimate's equivalence tests: Paths from sequences, higher-order
# network counts, and the multi-order model's likelihood. Moved from local_testing_and_equivalence/helper-python-equiv.R.
import pathpy as pp
import numpy as np

def _paths_from_sequences(seqs):
    paths = pp.Paths()
    for s in seqs:
        if len(s) >= 2:
            paths.add_path(','.join(list(s)), frequency=1)
    return paths

def hon_edges_df(paths, order):
    hon = pp.HigherOrderNetwork(paths, k=int(order))
    rows = []
    for (u, v), attrs in hon.edges.items():
        w = attrs.get('weight', np.zeros(2))
        # w is [subpath_count, longest_path_count]; Nestimate counts every occurrence
        rows.append((str(u), str(v), float(w[0]) + float(w[1])))
    return rows

def mogen_cumulative_likelihood(paths, max_order):
    # Total log-likelihood of MultiOrderModel with given max_order.
    # Matches Nestimate's $log_likelihood[[k+1]] at k = max_order,
    # which uses hierarchical decomposition across orders 0..max_order.
    mog = pp.MultiOrderModel(paths, max_order=int(max_order))
    return float(mog.likelihood(paths, log=True))

def mogen_hon_transition_matrix(paths, order):
    # Return (nodes, count_matrix) for a HigherOrderNetwork at given order.
    # Subpath + longest-path counts summed per cell. Node labels are pathpy's
    # comma-joined k-tuples.
    hon = pp.HigherOrderNetwork(paths, k=int(order))
    nodes = sorted(hon.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    M = np.zeros((n, n), dtype=float)
    for (u, v), attrs in hon.edges.items():
        w = attrs.get('weight', np.zeros(2))
        M[idx[u], idx[v]] = float(w[0]) + float(w[1])
    return {'nodes': nodes, 'counts': M.tolist()}

def mogen_stats(paths, max_order):
    mog = pp.MultiOrderModel(paths, max_order=int(max_order))
    # per-layer likelihoods, dof, and optimal order via LR test.
    # layer_likelihood(log=True) returns (log_lik, weight) tuple in pathpy 2.2.
    layers = []
    for k in range(int(max_order) + 1):
        lk_tuple = mog.layer_likelihood(paths, l=k, log=True)
        lk = float(lk_tuple[0]) if isinstance(lk_tuple, tuple) else float(lk_tuple)
        dof = mog.degrees_of_freedom(k)
        layers.append((int(k), lk, int(dof)))
    try:
        opt = int(mog.estimate_order(paths))
    except Exception:
        opt = -1
    return {'layers': layers, 'optimal_order': opt,
            'total_likelihood': float(mog.likelihood(paths, log=True)),
            'total_dof': int(mog.degrees_of_freedom(mog.max_order))}
