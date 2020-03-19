
'''

Draw protein domain-style plots, but containing the locations of the TEs

THis is the peptide-modified version;

The plots should be the mRNA, with a protein coding exon (if present), and the location of the TE;

This script collectes the data

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist, genome_sql
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.cm as cmap
from matplotlib.collections import PatchCollection

sys.path.append('../../../')
sys.path.append('../../')
import shared

draw = 'png'

def draw_domain(gene, filename, dfam, peptides):
    assert isinstance(peptides, list), 'peptides should be a list of tuples in the form [(l, r, peptide_seq), ...]'
    #print(gene)
    try:
        print('Doing %s' % gene['name'])
    except KeyError:
        print('Doing %s' % gene['enst'])

    if 'cds_local_locs' in gene:
        cdsl = gene['cds_local_locs'][0]
        cdsr = gene['cds_local_locs'][1]
    else:
        cdsl = None

    if cdsl == -1:
        cdsl = None

    _, tlength, splice_sites = shared.convert_genocode_to_splice_sites(gene)
    splice_sites = splice_sites[:-1] # last one is the termination;

    ts = gene['loc']['left']
    te = gene['loc']['right']

    fig = plot.figure(figsize=[3,1.6])
    fig.subplots_adjust(left=0.2, right=0.96, bottom=0.18, top=0.8)
    ax = fig.add_subplot(111)

    # three row markers
    line = []
    #line.append(mlines.Line2D([0, tlength], [1, 1], lw=1, alpha=1.0, color='grey', zorder=0))
    # TE lines:
    line.append(mlines.Line2D([0, tlength], [0.45, 0.45], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [0.15, 0.15], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [-0.15, -0.15], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [-0.45, -0.45], lw=1, alpha=1.0, color='grey', zorder=0))
    # gene line:
    line.append(mlines.Line2D([0, tlength], [-0.8, -0.8], lw=5, alpha=1.0, color='black', zorder=0))
    # Peptide line
    line.append(mlines.Line2D([0, tlength], [-1.1, -1.1], lw=1, alpha=1.0, color='grey', zorder=0))

    # TEs:
    patches = []
    this_dom = {}
    for dom_idx, dom in enumerate(gene['doms']):
        # get the type and subtype out of dfam:
        t = dfam.get(key='name', value=dom['dom'])[0]
        col_name = '%s:%s' % (t['type'], t['subtype'])
        color = shared.get_col(col_name)

        if t['type'] == 'LINE':
            posy = 0.45
        elif t['type'] == 'SINE':
            posy = 0.15
        elif t['type'] == 'LTR':
            posy=-0.15
        else: # And everything else?
            posy=-0.45

        s = dom['span'][0]
        e = dom['span'][1]
        if (s, e) in this_dom:
            continue # Just do one of the domains for clarity
        this_dom[(s, e)] = 1

        if dom['strand'] == '+':
            posx = s
            ha = 'left'
            strand = 0.06
        elif dom['strand'] == '-':
            posx = e
            ha = 'right'
            strand = -0.06

        #print(dom)
        # I need a color key for the TEs:
        rect = mpatches.Rectangle([s, posy], e-s, strand, ec="none", facecolor=color)
        ax.add_patch(rect)

        ax.text(posx, posy + (strand*2.1), dom['dom'], fontsize=5, ha=ha, va='center')

    # Peptides
    for p in peptides:
        s = p[0]
        e = p[1]
        if p[3] != 'No': # From a TE
            rect = mpatches.Rectangle([s, -1.1], e-s, 1.8, ec="none", alpha=0.5, facecolor='red')
        else:
            rect = mpatches.Rectangle([s, -1.1], e-s, 1.8, ec="none", alpha=0.5, facecolor='green')
        ax.add_patch(rect)

    # Draw splice markers;
    pad = (tlength / 120)
    for s in splice_sites:
        #line.append(mlines.Line2D([s-pad, s, s+pad], [-1.1, -1.0, -1.1], lw=1., alpha=0.8, color='r'))
        line.append(mlines.Line2D([s-pad, s+pad], [-0.85, -0.75], lw=1., alpha=0.8, color='r'))


    #del gene['doms']
    #print(gene)

    # CDS line:
    if cdsl is not None:
        cds_patch = mpatches.Rectangle([cdsl, -0.95], cdsr-cdsl, 0.3, ec="none", facecolor='black')
        ax.add_patch(cds_patch)

    # Finish drawing;
    collection = PatchCollection(patches)
    ax.add_collection(collection)
    [ax.add_line(l) for l in line]

    ax.set_xlim([0, tlength])
    ax.set_ylim([-1.2, 0.55])
    ax.set_yticks([-1.1, -0.8, -0.45, -0.15, 0, 0.15, 0.45])
    ax.set_yticklabels(['MS Peptides', 'Coding Seq', 'Retroposons', 'LTRs', 'TEs              ', 'SINEs', 'LINEs'])
    ax.set_title('%s(%s) %s %s' % (gene['enst'], gene['strand'], gene['transcript_id'], gene['name']), fontsize=6)

    ax.xaxis.set_tick_params(labelsize=6)
    ax.yaxis.set_tick_params(labelsize=6)

    ax.set_frame_on(False)
    ax.tick_params(left=False, bottom=False)

    if cdsl is not None:
        ax.set_xticks([0, cdsl, cdsr, tlength])
    else:
        ax.set_xticks([0, tlength])

    fig.savefig(filename)
    plot.close()

