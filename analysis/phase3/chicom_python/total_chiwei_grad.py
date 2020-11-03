from chiwei_grad import chiwei_grad
from math import sqrt
def total_chiwei_grad(iflag,p, weights_1, weights_sq_1, N_1, is_normalized_1, weights_2, weights_sq_2, N_2, is_normalized_2,grad_sum_old, s_old):
    grad_1 = chiwei_grad(p, weights_1, weights_sq_1, N_1, is_normalized_1)
    grad_2 = chiwei_grad(p, weights_2, weights_sq_2, N_2, is_normalized_2)
    if(grad_1==None or grad_2==None):
       return None, None, None
    bins = len(p)
    sumde=0.
    sumno=0.
    grad_total=[0.]*bins
    grad_sum = [grad_1[i] + grad_2[i] for i in range(bins)]
    coeff=0.
    if(iflag!=0):
       sumde= sum((grad_sum_old[i]*grad_sum_old[i] for i in range(bins)))  
       sumno= sum((grad_sum[i]*grad_sum[i]-grad_sum_old[i] for i in range(bins)))  
       if(sumde!=0.):
         coeff=sumno/sumde
       if((coeff<=0.) or (sumde==0.)):
         coeff=0.
    ggrad_sum_old = [grad_sum[i] for i in range(bins)] 
    grad_sum = [-grad_sum[i]+s_old[i]*coeff for i in range(bins)]
    ss_old = [grad_sum[i] for i in range(bins)]
    pp=sum(grad_sum[i] for i in range(bins))
    grad_total=[grad_sum[i]*bins-pp for i in range(bins)]
    normalize = sqrt(sum(grad_total[i]*grad_total[i] for i in range(bins)))
    grad_total = [grad_total[i] / normalize for i in range(bins)]
    return grad_total, ggrad_sum_old, ss_old

