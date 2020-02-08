import math

def y_heuristic_error(y_heuristic, rate_avg, rate_err, I_avg, I_err, P_avg, P_err):
    return y_heuristic*math.sqrt((rate_err/rate_avg)**2 + (I_err/I_avg)**2 + (P_err/P_avg)**2)

def y_heuristic(rate, I, P, Z):
    return rate/(I*P*Z**2)
