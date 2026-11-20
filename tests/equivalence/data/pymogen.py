# Reimplementation of MOGen model selection (Gote and Scholtes 2023: per-order log-likelihoods, degrees of
# freedom, AIC, BIC, likelihood-ratio tests) written for Nestimate's equivalence tests. Moved verbatim from
# local_testing_and_equivalence/test-mogen-pymogen-equiv.R. Standard library and scipy.stats.
import math
import sys
from collections import Counter, defaultdict
from scipy.stats import chi2

def build_mogen_ref(trajectories, max_order):
    '''Reference MOGen following Scholtes 2017 / Nestimate R implementation.

    Computes hierarchical log-likelihoods, DOF, AIC, BIC, LRT, and
    selects optimal order.
    '''

    # --- Step 0: unique first-order states ---
    all_states = sorted(set(s for traj in trajectories for s in traj))
    n_states = len(all_states)

    # --- Step 1: order-0 marginal distribution ---
    state_counts = Counter(s for traj in trajectories for s in traj)
    total = sum(state_counts.values())
    marginal = {s: state_counts[s] / total for s in all_states}

    # --- Step 2: build k-gram transition matrices for orders 1..max_order ---
    # Returns dict-of-dict: trans[src_tuple][tgt_tuple] = probability
    # Also count matrices for DOF computation
    def count_kgrams(trajs, k):
        '''Count transitions between consecutive k-grams.'''
        counts = defaultdict(Counter)
        for traj in trajs:
            n = len(traj)
            if n < k:
                continue
            # k-grams
            kgrams = [tuple(traj[i:i+k]) for i in range(n - k + 1)]
            # transitions between consecutive k-grams
            for i in range(len(kgrams) - 1):
                counts[kgrams[i]][kgrams[i+1]] += 1
        return counts

    def normalize(counts):
        '''Row-normalize count dict to transition probabilities.'''
        trans = {}
        for src, targets in counts.items():
            row_total = sum(targets.values())
            if row_total > 0:
                trans[src] = {t: c / row_total for t, c in targets.items()}
        return trans

    trans_mats = [marginal]  # index 0 = order 0 (marginal)
    count_mats = [None]       # placeholder for order 0
    for k in range(1, max_order + 1):
        cm = count_kgrams(trajectories, k)
        tm = normalize(cm)
        trans_mats.append(tm)
        count_mats.append(cm)

    # --- Step 3: DOF computation matching R's .mogen_layer_dof() ---
    # Order 0: n_states - 1
    # Order k >= 1: for each source context (row), count non-zero targets,
    #   DOF = max(n_nonzero - 1, 0), sum over all rows
    # But we need to count based on the TRANSITION MATRIX, not the raw counts.
    # In R, rows with zero row-sum have 0 DOF. Rows with row-sum > 0 are
    # normalized, and DOF = number of non-zero entries - 1.
    # Since we normalize (only keep rows with total > 0), every row in our
    # trans dict has at least 1 target. So DOF = sum(len(targets) - 1).
    # But R also includes rows that exist as nodes but have no outgoing edges
    # (row-sum = 0 => 0 DOF). Those rows don't appear in our dict, so they
    # contribute 0 automatically.
    layer_dofs = [n_states - 1]  # order 0
    for k in range(1, max_order + 1):
        tm = trans_mats[k]
        dof_k = sum(max(len(targets) - 1, 0) for targets in tm.values())
        layer_dofs.append(dof_k)

    cum_dofs = []
    running = 0
    for d in layer_dofs:
        running += d
        cum_dofs.append(running)

    # --- Step 4: log-likelihood (matching R's .mogen_log_likelihood) ---
    LOG_EPS = math.log(sys.float_info.epsilon)

    def log_likelihood(trajs, k, trans_mats_up_to_k):
        '''Compute LL of model order k.

        For each trajectory of length n:
          Step 1: log(marginal[state])
          Step i (i >= 2): order_used = min(i-1, k)
            if order_used == 0: log(marginal[state])
            else: look up trans_mats[order_used][src_kgram][tgt_kgram]

        This matches R's .mogen_log_likelihood exactly.
        '''
        p0 = trans_mats_up_to_k[0]  # marginal dict
        ll = 0.0
        for traj in trajs:
            n = len(traj)
            if n == 0:
                continue

            # Step 1: initial state
            prob = p0.get(traj[0], 0)
            ll += math.log(prob) if prob > 0 else LOG_EPS

            # Steps 2..n (1-indexed step i corresponds to 0-indexed position i-1)
            for pos in range(1, n):
                step = pos + 1  # 1-indexed step number
                order_used = min(step - 1, k)

                if order_used == 0:
                    prob = p0.get(traj[pos], 0)
                    ll += math.log(prob) if prob > 0 else LOG_EPS
                else:
                    # src_key: traj[pos - order_used : pos] (length = order_used)
                    # tgt_key: traj[pos - order_used + 1 : pos + 1] (length = order_used)
                    src_key = tuple(traj[pos - order_used : pos])
                    tgt_key = tuple(traj[pos - order_used + 1 : pos + 1])

                    tm = trans_mats_up_to_k[order_used]
                    prob = tm.get(src_key, {}).get(tgt_key, 0)
                    ll += math.log(prob) if prob > 0 else LOG_EPS
        return ll

    logliks = []
    for k in range(0, max_order + 1):
        ll = log_likelihood(trajectories, k, trans_mats[:k+1])
        logliks.append(ll)

    # --- Step 5: AIC, BIC ---
    n_obs = sum(len(t) for t in trajectories)
    aics = [2 * cum_dofs[k] - 2 * logliks[k] for k in range(max_order + 1)]
    bics = [math.log(n_obs) * cum_dofs[k] - 2 * logliks[k]
            for k in range(max_order + 1)]

    # --- Step 6: LRT (matching R's sequential test) ---
    # R uses layer_dof (not cumulative diff) as df_diff
    lrt_results = []
    for k in range(1, max_order + 1):
        x = -2 * (logliks[k-1] - logliks[k])
        df_diff = layer_dofs[k]
        if df_diff > 0 and x > 0:
            p_value = 1 - chi2.cdf(x, df_diff)
        else:
            p_value = 1.0
        lrt_results.append({
            'order': k, 'x': float(x),
            'delta_dof': df_diff, 'p': float(p_value)
        })

    # --- Step 7: Optimal order (AIC, BIC, LRT) ---
    opt_aic = int(aics.index(min(aics)))
    opt_bic = int(bics.index(min(bics)))

    opt_lrt = 0
    for res in lrt_results:
        if res['p'] < 0.01:
            opt_lrt = res['order']
        else:
            break

    return {
        'logliks': logliks,
        'layer_dofs': layer_dofs,
        'cum_dofs': cum_dofs,
        'aics': aics,
        'bics': bics,
        'lrt': lrt_results,
        'optimal_aic': opt_aic,
        'optimal_bic': opt_bic,
        'optimal_lrt': opt_lrt,
        'n_states': n_states,
        'n_obs': n_obs,
        'n_paths': len(trajectories)
    }
