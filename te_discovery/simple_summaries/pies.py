import glob, sys, os, gzip, numpy
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist

def lab_func(pct, allvals):
    absolute = int(pct/100.*sum(allvals))
    return "{:.1f}%\n{:,}".format(pct, absolute)

def pie(filename, data, labels, title=''):
    fig = plot.figure(figsize=[1,1])
    ax = fig.add_subplot(111)

    wedges, texts, autotexts = ax.pie(data, labels=labels, autopct=lambda pct: lab_func(pct, data))#, colors=['tomato', 'deepskyblue'])
    plot.setp(autotexts, size=6)
    plot.setp(texts, size=6)

    ax.set_title(title, size=6)
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    print('Saved %s' % filename)

def split_bar(filename, data_dict, title=''):

    # get all of the classes:
    all_keys = [] # preserve order
    for k in data_dict:
        for kk in data_dict[k]:
            if kk not in all_keys:
                all_keys.append(kk)
    print('Found {0} keys'.format(all_keys))

    vals = {k: [] for k in all_keys}

    labs = []
    for k in data_dict:
        labs.append(k)
        for kk in all_keys:
            vals[kk].append(float(data_dict[k][kk]))
    print(vals)

    scaled = {k: [] for k in all_keys}
    sums = None
    for k in all_keys:
        if sums is None:
            sums = numpy.zeros(len(vals[k]))
        sums += vals[k]

    for k in all_keys:
        vals[k] = numpy.array(vals[k])
        scaled[k] = numpy.array(vals[k])
        scaled[k] /= sums
        scaled[k] *= 100

    plot_hei = (0.8) - (0.05*len(labs))

    fig = plot.figure(figsize=[4,3])
    fig.subplots_adjust(left=0.35, right=0.95, bottom=plot_hei,)
    ax = fig.add_subplot(111)

    ypos = numpy.arange(len(data_dict))

    # data_dict = {'bar_row': {'class': 0, class2': 0}}

    bots = numpy.zeros(len(labs))
    for k in vals:
        ax.barh(ypos, scaled[k], 0.7, label=k, left=bots)
        for y, v, s, b in zip(ypos, vals[k], scaled[k], bots):
            ax.text(b+(s//2), y, '{0:,.0f} ({1:.0f}%)'.format(v, s), ha='center', va='center')
        bots += scaled[k]

    ax.set_yticks(ypos)
    ax.set_yticklabels(labs)

    ax.set_xlim([-2, 102])
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(['0%', '50%', '100%'])
    ax.set_title(title, size=6)
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    print('Saved %s' % filename)
