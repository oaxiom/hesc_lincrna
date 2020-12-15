import numpy as np
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from sklearn import linear_model
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error, r2_score

def scat(filename, x, y, xlabel, ylabel, xlims, ylims, alpha=0.1, cols=None, hlines=[0], vlines=[0]):
    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    xs = np.arange(min(x), max(x))
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    predict_y = intercept + slope * xs
    plot.plot(xs, predict_y, ':', c='black')
    ax.set_title('R={:.3f}; p={:.1e}'.format(r_value, p_value))

    ax.scatter(x, y, c=cols, s=4, alpha=alpha,  ec='none')

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    for l in hlines:
        ax.axvline(l, ls=":", lw=1.0, color="grey")
    for l in vlines:
        ax.axhline(l, ls=":", lw=1.0, color="grey")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)

def hist(filename, x, y, xlabel, ylabel, ranges, hlines=[0], vlines=[0]):
    # Hist2d:
    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    counts,xbins,ybins,image = ax.hist2d(x, y, #cmin=1,
        vmax=200,
        bins=60,
        range=ranges)

    plot.close(fig)

    fig = plot.figure(figsize=[2.4,2])
    ax = fig.add_subplot(111)

    xs = np.arange(-0.6, 0.6)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    predict_y = intercept + slope * xs

    plot.plot(xs, predict_y, ':', c='black')

    ax.set_title('R={:.3f}'.format(r_value))

    h = ax.hist2d(x, y, cmin=1,
        vmax=10, cmap=cm.plasma,
        bins=200, range=ranges)

    #ax.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
    #    linewidths=0.3,
    #    vmin=0,
    #    vmax=200,
    #    levels = [0, 1, 5, 10, 25, 50, 75, 100, 200],#, 300, 400]
    #    )

    ax.set_xlim(ranges[0])
    ax.set_ylim(ranges[1])

    for l in hlines:
        ax.axvline(l, ls=":", lw=1.0, color="grey")
    for l in vlines:
        ax.axhline(l, ls=":", lw=1.0, color="grey")

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
    ax.set_title('R={:.3f}; p={:.1e}'.format(r_value, p_value))

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


def boxplots(filename, data, qs,
    title=None,
    xlims=None,
    sizer=0.022,
    vert_height=4,
    bot_pad=0.1):

    mmheat_hei = 0.1+(sizer*len(data))
    fig = plot.figure(figsize=[2.8,vert_height])
    fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=bot_pad)
    ax = fig.add_subplot(111)
    ax.tick_params(right=True)

    m = 0
    ax.axvline(0, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
    ax.axvline(-0.25, ls=":", lw=0.5, color="grey")
    ax.axvline(+0.25, ls=":", lw=0.5, color="grey")

    dats = list(data.values())
    r = ax.boxplot(dats,
        showfliers=False,
        whis=True,
        patch_artist=True,
        widths=0.5, vert=False)

    plot.setp(r['medians'], color='black', lw=2) # set nicer colours
    plot.setp(r['boxes'], color='black', lw=0.5)
    plot.setp(r['caps'], color="grey", lw=0.5)
    plot.setp(r['whiskers'], color="grey", lw=0.5)

    ax.set_yticks(np.arange(len(data.keys()))+1)
    ax.set_yticklabels(data.keys())

    gtm = '#FF8A87'
    ltm = '#92A7FF'

    xlim = ax.get_xlim()[1]
    if xlims:
        ax.set_xlim(xlims)
        xlim = xlims[1]

    draw_qs = True

    for i, k, p in zip(range(0, len(data)), data, r['boxes']):
        if m >= 0.25:
            p.set_facecolor(gtm)
        elif m <= 0.00:
            p.set_facecolor(ltm)
        else:
            p.set_facecolor('lightgrey')

    if title:
        ax.set_title(title, fontsize=6)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)
    plot.close(fig)
