from total_chiwei import total_chiwei
from total_chiwei_grad import total_chiwei_grad
from step_size import step_size


def chicom(weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2, max_iter=1e+4):
    eps = 1e-10
    eps_xim = 1e-5
    eps_alp = 1e-5
    ilast=0
    maxic=5
    bins = len(weights_1)
    grad_sum_old=[0.]*bins
    s_old=[0.]*bins 

    if any(weights_1[i] == 0 and weights_2[i] == 0 for i in range(bins)):
        # No events in a bin in both histograms
        return 1, None, None, None

    ratio_1 = get_ratio(bins, weights_1, weights_sq_1)
    ratio_2 = get_ratio(bins, weights_2, weights_sq_2)

    weights_1_sum = sum(weights_1)
    weights_2_sum = sum(weights_2)

    # Initial estimate
    probs = [(weights_1[i] * ratio_1[i] + weights_2[i] * ratio_2[i]) / (ratio_1[i] * weights_1_sum + ratio_2[i] * weights_2_sum) for i in range(bins)]

    if any(probs[i] < eps for i in range(bins)):
        # Low statistics in a bin
        return 2, None, None, None

    current_best = 1e+300
    total_iters = 0
    iflag=0
    
    probs_grad = None

    ndf_1 = None
    ndf_2 = None

    for it in range(int(max_iter)):
        probs_grad, grad_sum_old, s_old= total_chiwei_grad(iflag, probs, weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2, grad_sum_old, s_old) 
        iflag=1
        total_iters = it + 1
        probs_ratios = [(-probs[i] / probs_grad[i]) for i in range(bins)]
        alpha_min = max(probs_ratios[i] for i in range(bins) if probs_ratios[i] < 0)
        alpha_max = min(probs_ratios[i] for i in range(bins) if probs_ratios[i] > 0)
        alpha_new  = step_size(probs, probs_grad, alpha_min, alpha_max, weights_1, weights_sq_1, N_1,
                              is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2)
        if(alpha_new==None):
             # Bad initial estimate
            return 3, None, None, None             
        stats, ndf_1, ndf_2, res_1, res_2 = total_chiwei(probs, probs_grad, alpha_new, alpha_min, alpha_max,
                                                         weights_1, weights_sq_1, N_1, is_normalized_1, weights_2,
                                                         weights_sq_2, N_2, is_normalized_2)
        if stats is None:
            # Bad initial estimate
            return 3, None, None, None
        delta = abs(current_best - stats) / stats
        if stats <= current_best:
           current_best = stats
        probs = [probs[i] + alpha_new * probs_grad[i] for i in range(bins)]
        if delta <= eps_xim and abs(alpha_new <= eps_alp):
           ilast=ilast+1
        else:
           ilast=0
        if ilast==maxic: 
                break
    if is_normalized_1 + is_normalized_2 == 2:
        ndf = bins - 1
    elif is_normalized_1 + is_normalized_2 == 1:
        ndf = bins - 2
        if ndf_1 == bins - 1 or ndf_2 == bins - 1:
            ndf = bins - 1
    else:
        ndf = bins - 2
        if ndf_1 == bins - 1 and ndf_2 == bins - 1:
            ndf = bins - 1
    stats=stats*bins 
    if(total_iters>=max_iter):
       return 4, stats, ndf, total_iters
    return 0, stats, ndf, total_iters


def get_ratio(bins, weights, weights_sq):
    ratio = [weights[i] / weights_sq[i] if weights[i] != 0 else 0 for i in range(bins)]
    default_ratio = sum(ratio) / sum(1 for r in ratio if r != 0)
    ratio = [default_ratio if r == 0 else r for r in ratio]
    return ratio
