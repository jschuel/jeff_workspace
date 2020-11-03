from math import sqrt


def chiwei(p, weights, weights_sq, n, is_normalized=False):
    """
    Calculates goodness of fit test statistics for weighted histograms.
    :param p: The probabilities for each bin. (list)
    :param weights: Total weight of events for each bin. (list)
    :param weights_sq: Total squared weight of events for each bin. (list)
    :param n: Number of events (integer)
    :param is_normalized: True if the histogram is normalized. Otherwise False. (Boolean)
    :return: The four following values. If error > 0, then all other values are None.
        error - Error code, 0 if the calculations are successful.
        stat - Test statistic,
        ndf - Number of degrees of freedom,
        residuals - True if calculations failed.
    """
    eps = 0.000001

    if any(value < 0 for value in p):
        # Data contains negative probability.
        return 101, None, None, None

    bins = len(p)

    if len(weights) is not bins:
        # The list 'weights' does not contain the same amount of weights as the number of bins.
        return 106, None, None, None

    if len(weights_sq) is not bins:
        # The list 'weights_sq' does not contain the same amount of weights as the number of bins.
        return 107, None, None, None

    p = [value / sum(p) for value in p]

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

    for i in range(bins):
        k = p_ratio.index(min(r for r in p_ratio if r is not None))
        p_ratio[k] = None

        sum_rwwp_k = sum_rwwp - ratio[i] * weights[i] * weights[i] / p[i]
        sum_rp_k = sum_rp - ratio[i] * p[i]
        sum_rw_k = sum_rw - ratio[i] * weights[i]

        c_k = 1
        if not is_normalized:
            c_k = sum_rp_k + sqrt(sum_rp_k / sum_rwwp_k) * (n - sum_rw_k)

        expected_frequencies = [n * p[j] * ratio[j] / c_k for j in range(bins)]
        if any(ef < 1 for ef in expected_frequencies):
            # Low statistics in a bin <1
            return 103, None, None, None

        if sum(1 for ef in expected_frequencies if ef < 5) > 0.2*bins:
            # Over 20% of bins with low statistics (<5)
            return 104, None, None, None

        if 1 - sum_rp_k / c_k >= eps:
            stat = c_k * sum_rwwp_k / n + pow(n - sum_rw_k, 2)/(n*(1 - sum_rp_k/c_k)) - n
            if stat >= 0:
                ndf = bins - 1 - (not is_normalized)
                res = [(weights[j] * c_k - n * p[j]) / sqrt(n * p[j] * (c_k / ratio[j] - p[j])) for j in range(bins)]
                return 0, stat, ndf, res

    # Correlation matrix non-positive
    return 102, None, None, None
