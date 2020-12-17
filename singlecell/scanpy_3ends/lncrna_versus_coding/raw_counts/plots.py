import os, sys, gzip, pickle
import matplotlib.pyplot as plot
import matplotlib.cm as cm
import numpy as np
from scipy.stats import linregress
from glbase3 import genelist, glload

def hist(filename, x, y, xlabel, ylabel, ranges, hlines=[0], vlines=[0]):
    # Hist2d:
    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    counts,xbins,ybins,image = ax.hist2d(x, y, #cmin=1,
        vmax=200,
        bins=60,
        range=ranges)

    plot.close(fig)

    fig = plot.figure(figsize=[4.4,4])
    spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[ 0.3, 1.0],
            height_ratios=[1.0, 0.3])

    ax = fig.add_subplot(spec[0,1])

    xs = np.arange(-0.6, 0.6)

    h = ax.hist2d(x, y,
        cmin=1,
        vmax=50,
        cmap=cm.plasma,
        bins=100,
        range=ranges
        )
    '''
    ax.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
        linewidths=0.3,
        vmin=0,
        vmax=200,
        levels = [0, 1, 5, 10, 25, 50, 75, 100, 200],#, 300, 400]
        )
    '''

    if ranges:
        ax.set_xlim(ranges[0])
        ax.set_ylim(ranges[1])

    for l in hlines:
        ax.axvline(l, ls=":", lw=1.0, color="grey")
    for l in vlines:
        ax.axhline(l, ls=":", lw=1.0, color="grey")

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    ax4 = fig.add_subplot(spec[1,0])
    plot.colorbar(h[3], ax=ax4)
    ax4.tick_params(left=None, bottom=None)
    ax4.set_xticklabels('')
    ax4.set_yticklabels('')

    ax = fig.add_subplot(spec[0, 0])
    ax.violinplot(ydata, [0], points=100, showmeans=True,)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    ax.set_ylabel(ylabel)
    ax.set_xticklabels('')
    ax.set_xticks([0])
    m = np.mean(ydata)
    ax.text(0, m, '{:.2f}'.format(m), fontsize=6)
    ax.set_ylim(ranges[1])

    ax = fig.add_subplot(spec[1,1])
    ax.violinplot(xdata, [0], points=100, vert=False, showmeans=True,)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    ax.set_xlabel(xlabel)
    ax.set_yticklabels('')
    ax.set_yticks([0])
    m = np.mean(xdata)
    ax.text(m, 0, '{:.2f}%'.format(m), fontsize=6)
    ax.set_xlim(ranges[0])

    fig.savefig(filename)


'''
res = {
    'm': {'coding': [], 'noncoding': []},
    's': {'coding': [], 'noncoding': []},
    'cv': {'coding': [], 'noncoding': []},
    'num_cells': {'coding': [], 'noncoding': []},
    }

'''

oh = open('res.pickle', 'rb')
res = pickle.load(oh)
oh.close()

for t in ['coding', 'noncoding']:

    xdata = res['num_cells'][t]
    ydata = np.log2(np.array(res['m'][t]))
    #ydata = res['s'][t]

    hist('hist-{}.pdf'.format(t), xdata, ydata, 'Percent of Cells', 'Mean of normalised expression',
        [[0, 50], [-0.20,8.1]], [], []
        )

