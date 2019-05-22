
'''

Draw protein domain-style plots, but containing the locations of the TEs

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

def draw_domain(gene, filename, gencode_db, dfam):
    #print(gene)
    try:
        print('Doing %s' % gene['name'])
    except KeyError:
        print('Doing %s' % gene['enst'])

    # Get cdsl if it has;
    gencode = None
    if 'ENST' in gene['enst']:
        gencode = gencode_db.find('%s (%s)' % (gene['name'].split(' ')[0], gene['enst']))
        #print(gencode, '%s (%s)' % (gene['name'].split(' ')[0], gene['enst']))
        if gencode:
            gencode = gencode[0]  # Should be 1 or 0 found
            # check it is not a non-coding one:
            if len(gencode['cds_loc']) == 0: # non-coding
                gencode = None

    #TODO: replace with shared.convert_genocode_to_local()
    cdsl = None; cdsr = None
    if gencode:
        # This is wrong if the transcrip is ~
        cdsl = gencode['cds_loc']['left']
        cdsr = gencode['cds_loc']['right']
        #print(gene['enst'], gene['transcript_id'], gene['name'], gene['strand'], gencode['loc'], gencode['cds_loc'], gene['exonStarts'], gene['exonEnds'])
    # Work out the mRNA length from the gene length and the exons and Es:
    tlength = 0
    currpos = gene['loc']['left'] # transcript position in genomic coords
    splice_sites = []
    newcdsl = 0 ; newcdsr = 0
    for splice in zip(gene['exonStarts'], gene['exonEnds']):
        if gencode:
            if cdsl >= splice[0] and cdsl <= splice[1]:
                newcdsl = tlength + (cdsl-splice[0])
            if cdsr >= splice[0] and cdsr <= splice[1]:
                newcdsr = tlength + (cdsr-splice[0])

        tlength += (splice[1]-splice[0])
        currpos = splice[1]
        splice_sites.append(tlength)

    if gencode and cdsl != cdsr: # if cdsl= cdsr then it is non-coding;
        cdsl = newcdsl
        cdsr = newcdsr

    if gene['strand'] == '-': # splice sites and cds are always + strand, so would need to invert them;
        splice_sites = [tlength-s for s in splice_sites]
        if gencode: # Have valid cdsl, cdsr
            cdsl = tlength - cdsl
            cdsr = tlength - cdsr

    splice_sites = splice_sites[:-1] # last one is the termination;

    ts = gene['loc']['left']
    te = gene['loc']['right']

    fig = plot.figure(figsize=[8,4])
    fig.subplots_adjust(left=0.2, right=0.95, bottom=0.18, top=0.8)
    ax = fig.add_subplot(111)

    # three row markers
    line = []
    line.append(mlines.Line2D([0, tlength], [1, 1], lw=1, alpha=1.0, color='grey', zorder=0))
    # TE lines:
    line.append(mlines.Line2D([0, tlength], [0.45, 0.45], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [0.15, 0.15], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [-0.15, -0.15], lw=1, alpha=1.0, color='grey', zorder=0))
    line.append(mlines.Line2D([0, tlength], [-0.45, -0.45], lw=1, alpha=1.0, color='grey', zorder=0))

    line.append(mlines.Line2D([0, tlength], [-1, -1], lw=5, alpha=1.0, color='black', zorder=0))

    # TEs:
    patches = []
    this_dom = {}
    for dom_idx, dom in enumerate(gene['doms']):
        # get the type and subtype out of dfam:
        t = dfam.get(key='name', value=dom['dom'])[0]
        col_name = '%s:%s' % (t['type'], t['subtype'])
        color = shared.col_keys[col_name]

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

        ax.text(posx, posy + (strand*2.1), dom['dom'], fontsize=6, ha=ha, va='center')

    # Draw splice markers;
    pad = (tlength / 120)
    for s in splice_sites:
        line.append(mlines.Line2D([s-pad, s, s+pad], [1.1, 1.0, 1.1], lw=1.5, alpha=0.8, color='r'))

    #del gene['doms']
    #print(gene)

    # CDS line:
    if cdsl:
        cds_patch = mpatches.Rectangle([cdsl, -1.2], cdsr-cdsl, 0.4, ec="none", facecolor='black')
        ax.add_patch(cds_patch)

    # Finish drawing;
    collection = PatchCollection(patches)
    ax.add_collection(collection)
    [ax.add_line(l) for l in line]

    ax.set_xlim([0, tlength])
    ax.set_ylim([-1.5, 1.5])
    ax.set_yticks([-1.0, -0.45, -0.15, 0, 0.15, 0.45, 1.0])
    ax.set_yticklabels(['Coding Seq', 'LINE', 'SINE', 'TEs              ', 'LTR', 'Other', 'Splice Junctions'])
    ax.set_title('%s(%s) %s %s' % (gene['enst'], gene['strand'], gene['transcript_id'], gene['name']))
    ax.set_frame_on(False)
    ax.tick_params(left=False)

    # set xticks:
    if gencode:
        ax.set_xticks([0, cdsl, cdsr, tlength])
    else:
        ax.set_xticks([0, tlength])


    fig.savefig(filename)
    plot.close()

