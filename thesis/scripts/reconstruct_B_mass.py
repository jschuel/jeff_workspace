import root_pandas as rp
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy

def get_data(input_file):
    data = rp.read_root(input_file, key='B')
    return data

def get_signal(data):
    sig = data.loc[data['isSignal'] == 1]
    return sig

def get_background(data):
    bg = data.loc[data['isSignal'] != 1]
    return bg


input_file = '~/workspace/thesis/signal_MC/signal_mc_exp55.root'
data = get_data(input_file)
sig = get_signal(data)
bg = get_background(data)
signal_to_background_fraction = data.isSignal.mean()
print(signal_to_background_fraction)
def pdf(B_mass, B_width, M):
    signal_pdf = scipy.stats.norm.pdf(M, loc=B_mass, scale=B_width)
    bckgrd_pdf = 1 #assume uniform background
    return signal_to_background_fraction * signal_pdf + (1-signal_to_background_fraction) * bckgrd_pdf

def logLikelihood(parameters, M):
    B_mass, B_width = parameters
    return -np.log(pdf(B_mass, B_width, M)).sum()

result = scipy.optimize.minimize(logLikelihood, [data.Mbc.mean(), data.Mbc.std()], args=(data.Mbc), method='Nelder-Mead')
print(f"Fit success: {result.success}, B mass: {result.x[0]}, width: {result.x[1]}")

plt.hist(sig['Mbc'], bins = 100, alpha = 1, histtype = 'stepfilled', lw = 1, range = (5.2, 5.3), edgecolor='black', label = 'signal')
plt.hist(bg['Mbc'], bins = 10, histtype = 'step', alpha = 1, range = (5.2, 5.3), label = 'background')
x = np.linspace(5.2, 5.3, 1000)
y = pdf(result['x'][0], result['x'][1], x) * len(data) * 0.5/450
plt.plot(x, y, color='tab:red', lw=1.5, label = 'fit: Mbc = %s +/- %s'%(float('%.4g' % result['x'][0]), float('%.1g' % result['x'][1])))
plt.legend()
plt.xlabel('Mbc [GeV]')
plt.grid()
plt.show()
