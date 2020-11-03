from math import sqrt


def chiwei_grad(p, weights, weights_sq, n, is_normalized):
    eps = 0.000001

    if any(value < 0 for value in p):
        # Data contains negative probability.
        return None

    bins = len(p)

    if len(weights) is not bins:
        # The list 'weights' does not contain the same amount of weights as the number of bins.
        return None

    if len(weights_sq) is not bins:
        # The list 'weights_sq' does not contain the same amount of weights as the number of bins.
        return None
    sum_p=sum(p) 
    p = [value / sum_p for value in p]

    ratio = [weights[i] / weights_sq[i] if weights[i] != 0 else 0 for i in range(bins)]
    default_ratio = sum(ratio) / sum(1 for r in ratio if r != 0)
    ratio = [default_ratio if r == 0 else r for r in ratio]

    # The algorithm requires calculating the sum for all i != k. Here the total sums are calculated,
    # and then the k'th element is removed in each iteration.
    sum_rwwp = sum([ratio[i] * weights[i] * weights[i] / p[i] for i in range(bins)])
    sum_rw = sum([ratio[i] * weights[i] for i in range(bins)])
    sum_rp = sum([ratio[i] * p[i] for i in range(bins)])

    # The indices are used in increasing order according to the p_ratio until success.
    p_ratio = [p[i] / ratio[i] for i in range(bins)]

    k = p_ratio.index(min(r for r in p_ratio if r is not None))
    p_ratio[k] = None

    sum_rwwp_k = sum_rwwp - ratio[k] * weights[k] * weights[k] / p[k]
    sum_rp_k = sum_rp - ratio[k] * p[k]
    sum_rw_k = sum_rw - ratio[k] * weights[k]

    c_k = 1
    if not is_normalized:
       c_k = sum_rp_k + sqrt(sum_rp_k / sum_rwwp_k) * (n - sum_rw_k)

    expected_frequencies = [n * p[j] * ratio[j] / c_k for j in range(bins)]
    if any(ef < 1 for ef in expected_frequencies):
        # Low statistics in a bin <1
       return None

    if sum(1 for ef in expected_frequencies if ef < 5) > 0.2*bins:
        # Over 20% of bins with low statistics (<5)
       return None

    if 1 - sum_rp_k / c_k >= eps:
        part_1 = [c_k * weights[i] * weights[i] / (p[i] * p[i]) for i in range(bins)]
        part_2 = (n - sum_rw_k) * (n - sum_rw_k) / (c_k * (1 - sum_rp_k / c_k) * (1 - sum_rp_k / c_k))
        grad = [ratio[i] * (part_2 - part_1[i]) / n if i != k else 0 for i in range(bins)]
        return grad

    # Correlation matrix non-positive
    return None
