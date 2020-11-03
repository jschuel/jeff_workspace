from total_chiwei import total_chiwei
from scipy.optimize import brent
from scipy.optimize import bracket
from scipy.optimize import minimize_scalar


def step_size(probs, probs_grad, alpha_min, alpha_max, weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2):

    f = lambda u: total_chiwei(probs, probs_grad, u, alpha_min, alpha_max, weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2)[0]
    toler=1.0e-05
    uu=abs(alpha_max)
    hh=abs(alpha_min)
    dd=min(uu,hh)
    h0=dd*toler 
    ffap=f(0.)
    ffbp=f(h0)
    if(ffap==None or ffbp==None):
       return None 
    xxa,xxb,xxc,fa,fb,fc,funcalls=bracket(f,xa=0.,xb=h0)
    if(fa==None or fb==None or fc==None):
       return None  
    if((fb<ffap) and abs(xxc-xxa)<=(toler*2.)):
       return xxb 
    res = minimize_scalar(f,bracket=(xxa, xxb, xxc), method='Golden')
    xmin=res.x
    if(xmin==None):
       return None 
    return xmin
