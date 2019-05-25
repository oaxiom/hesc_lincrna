import glob, sys, os, gzip
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
    fig.savefig(filename.replace('.png', '.svg'))
    print('Saved %s' % filename)
