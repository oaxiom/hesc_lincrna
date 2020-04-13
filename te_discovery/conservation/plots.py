import sys, os, itertools
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.tri as tri
from sklearn import linear_model
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error, r2_score
from glbase3 import *

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;

def scat(filename, x, y, xlabel, ylabel, xlims, ylims):
    fig = plot.figure(figsize=[3,3])
    ax = fig.add_subplot(111)

    xs = np.arange(-0.6, 0.6)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    predict_y = intercept + slope * xs

    plot.plot(xs, predict_y, ':', c='black')

    ax.set_title('R={:.3f}'.format(r_value))

    ax.scatter(x, y, s=4, alpha=0.1,  ec='none')

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    ax.axvline(0, ls=":", color="grey")
    ax.axhline(0, ls=":", color="grey")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)

def hist(filename, x, y, xlabel, ylabel, ranges):
    # Hist2d:
    fig = plot.figure(figsize=[3.4,3])
    ax = fig.add_subplot(111)

    xs = np.arange(-0.6, 0.6)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    predict_y = intercept + slope * xs

    plot.plot(xs, predict_y, ':', c='black')

    ax.set_title('R={:.3f}'.format(r_value))

    h = ax.hist2d(x, y, cmin=1,
        vmax=400,
        bins=20, range=ranges)

    ax.set_xlim(ranges[0])
    ax.set_ylim(ranges[1])

    ax.axvline(0, ls=":", color="grey")
    ax.axhline(0, ls=":", color="grey")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    plot.colorbar(h[3])

    fig.savefig(filename)

def contour(filename, x, y, xlabel, ylabel, ranges, vmax=200):
    # Hist2d:
    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    counts,xbins,ybins,image = ax.hist2d(x, y, #cmin=1,
        vmax=vmax,
        bins=60, range=ranges)

    plot.close(fig)

    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    ax.scatter(x, y, c='salmon', s=4, alpha=0.5,  ec='none')

    xs = np.arange(-0.6, 0.6)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    predict_y = intercept + slope * xs

    plot.plot(xs, predict_y, ':', c='black')

    ax.set_title('R={:.3f}'.format(r_value))

    ax.set_xlim(ranges[0])
    ax.set_ylim(ranges[1])

    ax.axvline(0, ls=":", color="grey")
    ax.axhline(0, ls=":", color="grey")

    ax.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
        linewidths=0.3,
        vmin=0, vmax=vmax,
        levels = [0, 1, 5, 10, 25, 50, 75, 100, 200],#, 300, 400]
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    plot.colorbar(image)
    fig.savefig(filename)

gl = glload('phyloP_conservation_table.glb')
print(gl)

scat('scat_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    xlims=[-0.6, 0.7],
    ylims=[-0.6, 0.7],
    )

hist('hist_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    ranges=[[-0.6, 0.7], [-0.6, 0.7]],
    )

contour('cont_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    ranges=[[-0.6, 0.7], [-0.6, 0.7]],
    )

for t in ('phyloP_tes', 'phyloP_nottes'):

    scat(filename='scat_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        xlims=[-0.6, 0.7],
        ylims=[-3, 9],
        )

    hist(filename='hist_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        ranges=[[-0.6, 0.7], [-3, 9]],
        )

    contour(filename='cont_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        ranges=[[-0.6, 0.7], [-3, 9]],
        vmax=100,
        )
