import numpy
import matplotlib.pyplot as plot
plot.rcParams['pdf.fonttype'] = 42

col_keys = {
    'DNA:.': 'black',
    'DNA:Crypton': 'royalblue',
    'DNA:Crypton-A': 'royalblue',
    'DNA:Kolobok': 'royalblue',
    'DNA:MULE-MuDR': 'skyblue',
    'DNA:Merlin': 'skyblue', ##
    'DNA:PIF-Harbinger': 'navy',
    'DNA:PiggyBac': 'navy',
    'DNA:TcMar': 'deepskyblue',
    'DNA:TcMar-Mariner': 'deepskyblue',
    'DNA:TcMar-Tc1': 'deepskyblue',
    'DNA:TcMar-Tc2': 'deepskyblue',
    'DNA:TcMar-Tigger': 'deepskyblue',

    'DNA:hAT': 'darkslategrey',
    'DNA:hAT-Ac': 'darkslategrey',
    'DNA:hAT-Blackjack': 'darkslategrey',
    'DNA:hAT-Charlie': 'darkslategrey',
    'DNA:hAT-Tag1': 'darkslategrey',
    'DNA:hAT-Tip100': 'darkslategrey',

    'LINE:CR1': 'khaki',
    'LINE:I-Jockey': 'gold',
    'LINE:L1': 'goldenrod',
    'LINE:L2': 'goldenrod',
    'LINE:Penelope': 'goldenrod',
    'LINE:RTE-BovB': 'goldenrod',

    'LTR:.': 'black',
    'LTR:ERV1': 'tomato',
    'LTR:ERVK': 'orangered',
    'LTR:ERVL': 'darkred',
    'LTR:ERVL-MaLR': 'brown',
    'LTR:Gypsy': 'indianred',
    'LINE:RTE-X': 'coral',

    'RC:Helitron': 'greenyellow',

    'Retroposon:SVA': 'orange', # To get new cols;

    'Satellite:.': 'black',
    'Satellite:acromeric': 'olive',
    'Satellite:centromeric': 'olivedrab',
    'Satellite:subtelomeric': 'olivedrab',

    'SINE:Alu': 'black',
    'SINE:5S-Deu-L2': 'black',
    'SINE:tRNA-Deu': 'black',
    'SINE:tRNA': 'black',
    'SINE:MIR': 'black',

    'Unknown:.': 'grey',

    'rRNA:.': 'purple',
    'scRNA:.': 'blueviolet',
    'snRNA:.': 'violet',
    'tRNA:.': 'orchid',
    }

def get_col(e):
    if e in col_keys:
        return col_keys[e]
    print('Colors: %s not found' % e)
    return 'grey'

def get_cols(labels):
    return [get_col(e) for e in labels]

part1 = {'~)': 'variant-isoform',
    '=)': 'known-isoform',
    '!)': 'unknown-isoform',}

part2 = {'(ME': 'multi-exon',
    '(SE': 'single-exon'}

part3 = {'C': 'coding',
    'NC': 'noncoding',
    'U': 'unknown'}

part4 = {'SR': 'short-read',
    'SR+LR': 'long-read', # Takes preference for purposes of proof
    'LR': 'long-read'}

part5 = {'ES+': 'ES-enriched',
    'ES:': 'ES-neutral',
    'ES-': 'ES-depleted'}

def classify_transcript(name):
    sum = name.split(' ')[1].split(';')
    if 'HSC' in name: # and hence ';!)'
        destination = 'novel_transcript_%s_%s_%s_%s' % (part3[sum[1]], part2[sum[0]], part4[sum[3]], part5[sum[2]])
        alpha = '.'

    else:
        # transcript_type - coding - exon_status - evidence - expression
        destination = '%s_%s_%s_%s_%s' % (part1[sum[4]], part3[sum[1]], part2[sum[0]], part4[sum[3]], part5[sum[2]])

        alpha = name[0]

    return(destination, alpha)

def convert_genocode_to_local(gencode):
    '''
    Convert a gencode genomic annotation into a local transcript structure. i.e. convert this:

    [{'loc': <location chr1:57598-64116>,
    'cds_loc': <location chr1:57598-57598>,
    'exonStarts': [57598, 58700, 62916],
    'exonEnds': [57653, 58856, 64116],
    'name': 'OR4G11P (ENST00000642116)',
    'type': 'gene',
    'strand': '+'}]

    into:
    TSS, TTS, CDSL, CDSR, splice_locations

    and it is always on the 5' strand, and TSS (by definition) == 0

    '''
    global_cdsl = gencode['cds_loc']['left']
    global_cdsr = gencode['cds_loc']['right']
    global_left = gencode['loc']['left']
    global_right = gencode['loc']['right']
    local_cdsl = -1
    local_cdsr = -1

    # These are not garunteed to be in the correct order in the input:
    exonStarts = sorted(gencode['exonStarts'])
    exonEnds = sorted(gencode['exonEnds'])
    exonCount = len(gencode['exonEnds'])

    splice_sites = []
    newcdsl = 0 ; newcdsr = 0

    # convert all of the positions into the spliced locations, on the + strand;
    if gencode['strand'] == '+':
        tlength = 0
        for exonidx, exon in enumerate(zip(exonStarts, exonEnds)):
            if global_cdsl >= exon[0] and global_cdsl <= exon[1]:
                local_cdsl = tlength + (global_cdsl-exon[0])

            if global_cdsr >= exon[0] and global_cdsr <= exon[1]:
                local_cdsr = (tlength+3+exonCount) + (global_cdsr-exon[0]) # strangeness for 0-based versus 1-based and open/closed :(

            #print(exon, global_cdsl, global_cdsr, local_cdsl, local_cdsr, tlength)
            tlength += (exon[1]-exon[0])
            splice_sites.append(tlength)

    elif gencode['strand'] == '-':
        tlength = 0
        exonStarts = reversed(exonStarts)
        exonEnds = reversed(exonEnds)
        for exon in zip(exonStarts, exonEnds):
            if global_cdsr >= exon[0] and global_cdsr <= exon[1]:
                local_cdsl = (exon[1]-global_cdsr) + tlength

            if global_cdsl >= exon[0] and global_cdsl <= exon[1]:
                local_cdsr = (exon[1]-global_cdsl) + (tlength+3+exonCount) # strangeness for 0-based versus 1-based and open/closed :(

            #print(exon, global_cdsl, global_cdsr, local_cdsl, local_cdsr, tlength)
            tlength += (exon[1]-exon[0])
            splice_sites.append(tlength)

        # check estimate is the same:
        if tlength != get_transcript_length(gencode):
            print(tlength,get_transcript_length(gencode))
            1/0
    tlength = get_transcript_length(gencode)

    splice_sites = splice_sites[:-1] # last one is the termination;

    #ts = gene['loc']['left']
    #te = gene['loc']['right']

    return 0, tlength, local_cdsl, local_cdsr, splice_sites

def convert_genocode_to_splice_sites(gencode):
    '''
    Convert a gencode genomic annotation into a local transcript structure. i.e. convert this:

    [{'loc': <location chr1:57598-64116>,
    'cds_loc': <location chr1:57598-57598>,
    'exonStarts': [57598, 58700, 62916],
    'exonEnds': [57653, 58856, 64116],
    'name': 'OR4G11P (ENST00000642116)',
    'type': 'gene',
    'strand': '+'}]

    into:
    TSS, TTS, splice_locations

    and it is always on the 5' strand, and TSS (by definition) == 0

    '''
    # Work out the mRNA length from the gene length and the exons and Es:
    tlength = 0
    currpos = gencode['loc']['left'] # transcript position in genomic coords
    splice_sites = []
    for splice in zip(gencode['exonStarts'], gencode['exonEnds']):
        tlength += (splice[1]-splice[0])
        currpos = splice[1]
        splice_sites.append(tlength)

    splice_sites = splice_sites[:-1] # last one is the termination;

    return 0, tlength, splice_sites

def get_transcript_length(gencode):
    '''
    Report the length of the transcript, from something that looks like this:

    [{'loc': <location chr1:57598-64116>,
    'exonStarts': [57598, 58700, 62916],
    'exonEnds': [57653, 58856, 64116],
    'name': 'OR4G11P (ENST00000642116)',
    'type': 'gene',
    'strand': '+'}]

    and it is always on the 5' strand, and TSS (by definition) == 0

    '''
    # Work out the mRNA length from the gene length and the exons and Es:
    tlength = 0
    for splice in zip(gencode['exonStarts'], gencode['exonEnds']):
        tlength += (splice[1]-splice[0])

    tlength = tlength + (gencode['loc']['right'] - splice[1]) # add the last segment

    return tlength

def pickle_it(filename, object):
    import pickle

    oh = open(filename, 'wb')
    pickle.dump(object, oh, protocol=4)
    oh.close()

def get_pickle(filename):
    import pickle

    oh = open(filename, 'rb')
    res = pickle.load(oh)
    oh.close()
    return res

table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

def split3(s):
    return [s[i:i+3] for i in range(0, len(s), 3)]

def translateAA(seq):
    return [table[codon.upper()] for codon in split3(seq)]

def nice_scatter(x=None, y=None, filename=None, do_best_fit_line=False, spot_cols='grey',
    print_correlation=False, spot_size=4, label_fontsize=14, label=False, label_t=1.301,
    label_tester=None,
    **kargs):
    """
    **Purpose**
        Draw a nice simple scatter plot

        Extended from glbase3
    """
    from glbase3 import draw

    draw = draw()

    fig = draw.getfigure(**kargs)
    ax = fig.add_subplot(111)

    ax.scatter(x, y, s=spot_size, c=spot_cols, alpha=0.2, edgecolors="none")

    if label:
        for x, y, t in zip(x, y, label):
            if label_tester:
                if True in [test in t for test in label_tester]:
                    ax.text(x, y, t, fontsize=5)
                continue
            if y > label_t: # q=0.05
                ax.text(x, y, t, fontsize=5)

    draw.do_common_args(ax, **kargs)

    return(draw.savefigure(fig, filename))

if __name__ == '__main__':
    import matplotlib.pyplot as plot

    # legend:
    fig, ax=plot.subplots(figsize=(4,6))
    ax.set_xlim(0, 30)
    ax.set_axis_off()

    for i, name in enumerate(reversed(list(col_keys.keys()))):
        y = i * 10

        ax.text(2.5, y, name, fontsize=6,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, 1, 2, color=col_keys[name], linewidth=5)
    fig.savefig('legend.png')
    fig.savefig('legend.svg')
    fig.savefig('legend.pdf')

    print('\n'.join(list(col_keys.keys())))

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


def lab_func(pct, allvals):
    absolute = int(pct/100.*sum(allvals))
    return "{:.1f}%\n{:,}".format(pct, absolute)

def split_bar(filename, data_dict, key_order=None, title='', cols=None, figsize=[4,3]):
    if not cols:
        cols = plot.rcParams['axes.prop_cycle'].by_key()['color']

    # get all of the classes:
    if not key_order:
        all_keys = [] # preserve order
        for k in data_dict:
            for kk in data_dict[k]:
                if kk not in all_keys:
                    all_keys.append(kk)
        print('Found {0} keys'.format(all_keys))
    else:
        all_keys = key_order

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

    plot_hei = (0.8) - (0.04*len(labs))

    plot.rcParams['pdf.fonttype'] = 42
    fig = plot.figure(figsize=[4,3])
    fig.subplots_adjust(left=0.35, right=0.95, bottom=plot_hei,)
    ax = fig.add_subplot(111)
    ax.set_prop_cycle('color', cols)

    ypos = numpy.arange(len(data_dict))

    # data_dict = {'bar_row': {'class': 0, class2': 0}}

    bots = numpy.zeros(len(labs))
    for k in vals:
        ax.barh(ypos, scaled[k], 0.7, label=k, left=bots)
        for y, v, s, b in zip(ypos, vals[k], scaled[k], bots):
            ax.text(b+(s//2), y, '{0:,.0f} ({1:.0f}%)'.format(v, s), ha='center', va='center', fontsize=6)
        bots += scaled[k]

    ax.set_yticks(ypos)
    ax.set_yticklabels(labs)

    ax.set_xlim([-2, 102])
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(['0%', '50%', '100%'])
    ax.set_title(title, size=6)
    ax.legend()
    plot.legend(loc='upper left', bbox_to_anchor=(0.0, -0.4), prop={'size': 6})
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    print('Saved %s' % filename)

def bar(filename, data_dict, key_order=None, title='', cols=None, figsize=[4,3]):
    if not cols:
        cols = plot.rcParams['axes.prop_cycle'].by_key()['color']

    # get all of the classes:
    if not key_order:
        all_keys = [] # preserve order
        for k in data_dict:
            for kk in data_dict[k]:
                if kk not in all_keys:
                    all_keys.append(kk)
        print('Found {0} keys'.format(all_keys))
    else:
        all_keys = key_order

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

    plot_hei = (0.8) - (0.04*len(labs))

    plot.rcParams['pdf.fonttype'] = 42
    fig = plot.figure(figsize=[4,3])
    fig.subplots_adjust(left=0.35, right=0.95, bottom=plot_hei,)
    ax = fig.add_subplot(111)
    ax.set_prop_cycle('color', cols)

    ypos = numpy.arange(len(data_dict))

    # data_dict = {'bar_row': {'class': 0, class2': 0}}

    bots = numpy.zeros(len(labs))
    for k in vals:
        ax.barh(ypos, vals[k], 0.7, label=k, left=bots)
        for y, v, s, b in zip(ypos, vals[k], vals[k], bots):
            ax.text(b+(s//2), y, '{0:,.0f} ({1:.0f}%)'.format(v, s), ha='center', va='center', fontsize=6)
        bots += vals[k]

    ax.set_yticks(ypos)
    ax.set_yticklabels(labs)

    ax.set_title(title, size=6)
    ax.legend()
    plot.legend(loc='upper left', bbox_to_anchor=(0.0, -0.4), prop={'size': 6})
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    print('Saved %s' % filename)

def boxplots(filename, data, qs, no_TE_key='no TE'):
    fig = plot.figure(figsize=[2.8,4])
    mmheat_hei = 0.1+(0.022*len(data))
    fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=0.1)
    ax = fig.add_subplot(111)
    ax.tick_params(right=True)

    if no_TE_key:
        m = numpy.median(data[no_TE_key])
        ax.axvline(m, ls=":", lw=0.5, color="grey") # add a grey line at m for better orientation
        ax.axvline(m-1, ls=":", lw=0.5, color="grey")
        ax.axvline(m+1, ls=":", lw=0.5, color="grey")
    else: # Probably the fold-change ones
        m = 0
        ax.axvline(0, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
        ax.axvline(-1, ls=":", lw=0.5, color="grey")
        ax.axvline(+1, ls=":", lw=0.5, color="grey")

    dats = numpy.array(list(data.values()))
    r = ax.boxplot(dats,
        showfliers=False,
        whis=True,
        patch_artist=True,
        widths=0.5, vert=False)

    plot.setp(r['medians'], color='black', lw=2) # set nicer colours
    plot.setp(r['boxes'], color='black', lw=0.5)
    plot.setp(r['caps'], color="grey", lw=0.5)
    plot.setp(r['whiskers'], color="grey", lw=0.5)

    ax.set_yticklabels(data.keys())

    gtm = '#FF8A87'
    ltm = '#92A7FF'

    xlim = ax.get_xlim()[1]

    draw_qs = True
    if not qs:
        # Change the qs to enforce some kind of criteria, like 2-fold.
        qs = {}
        for k in data:
            m = numpy.median(data[k])
            if abs(m) > 0.53:
                qs[k] = 0.00001# spoof to force colour drawing;
            else:
                qs[k] = 1  # i.e. grey
        draw_qs = False # But don't draw;

    for i, k, p in zip(range(0, len(data)), data, r['boxes']):
        #ax.text(6.3, i+0.5, q, ha='left', va='center', fontsize=6,)
        print(k, qs[k])
        if qs[k] < 0.05:
            if draw_qs: ax.text(xlim+(xlim/5), i+1, '*{0:.1e}'.format(qs[k]), ha='left', va='center', fontsize=6,)
            if numpy.median(data[k]) > m:
                p.set_facecolor(gtm)
            else:
                p.set_facecolor(ltm)
        else:
            if draw_qs: ax.text(xlim+(xlim/5), i+1, '{0:.1e}'.format(qs[k]), ha='left', va='center', fontsize=6,)
            p.set_facecolor('lightgrey')

        if k == 'no TE':
            p.set_facecolor('grey')

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)
