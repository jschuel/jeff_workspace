from chiwei import chiwei


def total_chiwei(p, p_grad, alpha, alpha_min, alpha_max, weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2):
    bins = len(weights_1) 
    if alpha > alpha_max or alpha < alpha_min:
        return None, None, None, None, None
    p_new = [p[i] + alpha * p_grad[i] for i in range(len(p))]
    error_1, stat_1, ndf_1, res_1 = chiwei(p_new, weights_1, weights_sq_1, N_1, is_normalized_1)
    error_2, stat_2, ndf_2, res_2 = chiwei(p_new, weights_2, weights_sq_2, N_2, is_normalized_2)
    chiwei_sum = None
    if error_1 == 0 and error_2 == 0:
        chiwei_sum = (stat_1 + stat_2)/bins

    return chiwei_sum, ndf_1, ndf_2, res_1, res_2
