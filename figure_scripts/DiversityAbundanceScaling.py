from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table


p, fr, _lw, w, fs, sz = 2, 0.75, 0.5, 1, 6, 2
smin = False

mydir = os.path.expanduser('~/GitHub/DeepLife_16S')
sys.path.append(mydir+'/DiversityTools')
mydir2 = os.path.expanduser("~/")



def figplot(x, y, xlab, ylab, fig, n):

    fig.add_subplot(2, 2, n)
    plt.xscale('log')
    plt.yscale('log')
    plt.minorticks_off()

    d = pd.DataFrame({'x': np.log10(x)})
    d['y'] = np.log10(y)
    f = smf.ols('y ~ x', d).fit()

    m, b, r, p, std_err = stats.linregress(np.log10(x), np.log10(y))
    st, data, ss2 = summary_table(f, alpha=0.05)
    fitted = data[:,2]
    mean_ci_low, mean_ci_upp = data[:,4:6].T
    ci_low, ci_upp = data[:,6:8].T

    x, y, fitted, ci_low, ci_upp = zip(*sorted(zip(x, y, fitted, ci_low, ci_upp)))

    x = np.array(x)
    y = np.array(y)
    fitted = 10**np.array(fitted)
    ci_low = 10**np.array(ci_low)
    ci_upp = 10**np.array(ci_upp)

    if n == 1: lbl = r'$rarity$'+ ' = '+str(round(10**b,1))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'
    elif n == 2: lbl = r'$Nmax$'+ ' = '+str(round(10**b,1))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'
    elif n == 3: lbl = r'$Ev$'+ ' = '+str(round(10**b,1))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'
    elif n == 4: lbl = r'$S$'+ ' = '+str(round(10**b,1))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'

    plt.scatter(x, y, s = sz, color='0.2', linewidths=0.0, edgecolor=None)
    plt.fill_between(x, ci_upp, ci_low, color='0.5', lw=0.1, alpha=0.2)
    plt.plot(x, fitted,  color='k', ls='--', lw=0.5, label = lbl)

    plt.legend(loc='best', fontsize=8, frameon=False)

    plt.xlabel(xlab, fontsize=10)
    plt.ylabel(ylab, fontsize=10)
    plt.tick_params(axis='both', labelsize=6)
    
    return fig

file_name = 'GPC_OTU_Table_97sim_Sample_Div_Measures.txt'
#file_name = 'Louca-SADMetricData.txt'

df = pd.read_csv(mydir + '/diversity_data/'+file_name, delimiter=r"\s+")
df = df[df['N'] > 0]

df2 = pd.DataFrame({'N' : df['N']})

df2['D'] = df['Nmax']
df2['S'] = df['S']
df2['E'] = df['ESimp']
df2['R'] = df['logmodskew'] + 0.5

if smin: df2 = df2[df2['S'] > 1]
df2 = df2.replace([np.inf, -np.inf], np.nan).dropna()


fig = plt.figure()

xlab = '$N$'
ylab = 'Rarity'
fig = figplot(df2['N'], df2['R'], xlab, ylab, fig, 1)

xlab = '$N$'
ylab = 'Dominance'
fig = figplot( df2['N'], df2['D'], xlab, ylab, fig, 2)

xlab = '$N$'
ylab = 'Evenness'
fig = figplot(df2['N'], df2['E'], xlab, ylab, fig, 3)

xlab = '$N$'
ylab = 'Richness'
fig = figplot(df2['N'], df2['S'], xlab, ylab, fig, 4)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/figures/DiversityAbundanceScaling.png', dpi=400, bbox_inches = "tight")
plt.close()