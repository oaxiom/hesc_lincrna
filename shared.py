import numpy
import matplotlib.pyplot as plot
from sklearn import linear_model
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error, r2_score
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

def convert_local_to_genome(local_left, local_rite, loc, exonStarts, exonEnds, strand):
    '''
    convert local_left and local_rite to genomic coordinates;


    '''

    genome_left = 0
    genome_rite = 0

    loc_left = loc.loc['left']
    loc_rite = loc.loc['right']

    intron_length = 0
    exon_length = 0
    last_exon_rite = loc_left

    #print(loc)
    # convert all of the positions into the spliced locations, on the + strand;

    if strand == '-': # I think this doesn't matter, if '-' then you can just swap the locs?
        exonStarts.reverse()
        exonEnds.reverse()

    for exonidx, exon in enumerate(zip(exonStarts, exonEnds)):
        intron_length += (exon[0] - last_exon_rite)
        last_exon_rite = exon[1]

        local_exon_left = exon[0]- loc_left - intron_length
        local_exon_rite = exon[1] - intron_length - loc_left

        if local_left >= local_exon_left and local_left <= local_exon_rite: # within this exon;
            genome_left = exon[0] + (local_left-local_exon_left)

        if local_rite >= local_exon_left and local_rite <= local_exon_rite: # within this exon;
            genome_rite = exon[0] + (local_rite-local_exon_left)

        print(exonidx, exon, 'tlength', exon_length, intron_length, last_exon_rite, 'local_exon', local_exon_left, local_exon_rite, genome_left, genome_rite)
        exon_length += (exon[1]-exon[0])


    if strand == '+':
        pass

    # check estimate is the same:
    if exon_length+intron_length != loc_rite-loc_left:
        print('Size mismatch!', exon_length+intron_length, loc_rite - loc_left)
        #1/0

    return (genome_left, min([genome_rite, last_exon_rite]))

if __name__ == '__main__':
    from glbase3 import location
    print(convert_local_to_genome(289, 3373, location(loc='chr10:21534232-21743630'), # ENST00000377072.8
        [21534232, 21534645, 21538833, 21586294, 21595331, 21612348, 21614831, 21617112, 21651673, 21670449, 21673350, 21681332, 21682225, 21688491, 21713772, 21726244, 21727856, 21730900, 21732899, 21733504, 21733768, 21735139, 21740154, 21741939],
        [21534520, 21534804, 21538912, 21586348, 21595440, 21612451, 21614924, 21617207, 21651768, 21670704, 21673919, 21681376, 21682257, 21688538, 21713950, 21726355, 21727928, 21731054, 21733087, 21733592, 21734129, 21735235, 21740236, 21743630],
        '+')) #
    print('Answer: chr10:21534645-21740231')
    print(convert_local_to_genome(419, 2126, location(loc='chr10:3104724-3137718'), # ENST00000415005.6
        [3104724, 3105393, 3107214, 3108701, 3109355, 3112222, 3113119, 3113372, 3116776, 3118782, 3119892, 3129819, 3132380, 3133203, 3134483, 3135736, 3136450],
        [3105159, 3105501, 3107309, 3108793, 3109480, 3112286, 3113188, 3113518, 3116846, 3118869, 3120044, 3129983, 3132441, 3133314, 3134582, 3135838, 3137718],
        '+'))
    print('Answer: chr10:3105143-3136576')

    print(convert_local_to_genome(420, 2706, location(loc='chr9:119166629-119369435'), # 'ENST00000265922.8'
        [119369056, 119313138, 119248960, 119242047, 119238655, 119213919, 119208719, 119166629],
        [119369435, 119313405, 119249150, 119242216, 119238760, 119214155, 119208941, 119168224],
        '-'))
    print('Answer: chr9:119167087-119313355')

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
    label_tester=None, doR=False,
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

    if doR:
        xs = numpy.arange(min(x), max(x))
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        predict_y = intercept + slope * xs
        plot.plot(xs, predict_y, ':', c='black')
        ax.set_title('R={:.3f}; p={:.1e}'.format(r_value, p_value))

    ax.scatter(x, y, s=spot_size, c=spot_cols, alpha=0.2, edgecolors="none")

    if label:
        for x, y, t in zip(x, y, label):
            if label_tester:
                if True in [test in t for test in label_tester]:
                    if t: ax.text(x, y, t, fontsize=5)
                continue
            if y > label_t: # q=0.05
                ax.text(x, y, t, fontsize=5)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    draw.do_common_args(ax, **kargs)

    return(draw.savefigure(fig, filename))

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
        print('Found {} keys'.format(all_keys))
    else:
        all_keys = key_order

    vals = {k: [] for k in all_keys}

    labs = []
    for k in data_dict:
        labs.append(k)
        for kk in all_keys:
            vals[kk].append(float(data_dict[k][kk]))
    print('vals:', vals)

    # data_dict = {'bar_row': {'class': 0, class2': 0}}
    percs = {}
    for k in vals:
        s = sum(vals[k])
        percs[k] = [((n / s) * 100.0) for n in vals[k]]

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

    bots = numpy.zeros(len(labs))
    for k in vals:
        ax.barh(ypos, vals[k], 0.7, label=k, left=bots)
        for y, v, s, b in zip(ypos, vals[k], percs[k], bots):
            ax.text(b+(v//2), y, '{:,.0f} ({:.0f}%)'.format(v, s), ha='center', va='center', fontsize=6)
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

def boxplots(filename, data, qs, no_TE_key='no TE', title=None,
    xlims=None,
    trim_low_samples=False,
    sizer=0.022,
    vert_height=4,
    bot_pad=0.1):

    mmheat_hei = 0.1+(sizer*len(data))
    fig = plot.figure(figsize=[2.8,vert_height])
    fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=bot_pad)
    ax = fig.add_subplot(111)
    ax.tick_params(right=True)

    if no_TE_key:
        m = numpy.median(data[no_TE_key])
        ax.axvline(m, ls=":", lw=0.5, color="grey") # add a grey line at m for better orientation
        ax.axvline(m-1, ls=":", lw=0.5, color="grey")
        ax.axvline(m+1, ls=":", lw=0.5, color="grey")
    elif qs == '100%quartiles':
        ax.axvline(25, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
        ax.axvline(50, ls=":", lw=0.5, color="grey")
        ax.axvline(75, ls=":", lw=0.5, color="grey")
    else: # Probably the fold-change ones
        m = 0
        ax.axvline(0, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
        ax.axvline(-1, ls=":", lw=0.5, color="grey")
        ax.axvline(+1, ls=":", lw=0.5, color="grey")

    # blank out really empty ones:
    if trim_low_samples:
        newd = {}
        for k in data:
            if len(data[k]) <= trim_low_samples:
                newd[k] = []
            else:
                newd[k] = data[k]
        data = newd

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

    ax.set_yticklabels(data.keys())

    gtm = '#FF8A87'
    ltm = '#92A7FF'

    xlim = ax.get_xlim()[1]
    if xlims:
        ax.set_xlim(xlims)
        xlim = xlims[1]

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

    if qs == '100%quartiles':
        for i, k, p in zip(range(0, len(data)), data, r['boxes']):
            m = numpy.median(data[k])
            if m >= 75:
                p.set_facecolor(gtm)
            elif m <= 25:
                p.set_facecolor(ltm)
            else:
                p.set_facecolor('lightgrey')

    else:
        for i, k, p in zip(range(0, len(data)), data, r['boxes']):
            #ax.text(6.3, i+0.5, q, ha='left', va='center', fontsize=6,)
            #print(k, qs[k])
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

    if title: ax.set_title(title, fontsize=6)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)
    plot.close(fig)

def boxplots_simple(filename, data, qs,
    title=None,
    xlims=None,
    sizer=0.022,
    vert_height=4,
    col='lightgrey',
    bot_pad=0.1,
    vlines=[]):

    mmheat_hei = 0.1+(sizer*len(data))
    fig = plot.figure(figsize=[2.8,vert_height])
    fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=bot_pad)
    ax = fig.add_subplot(111)
    ax.tick_params(right=True)

    m = 0
    for vli in vlines:
        ax.axvline(vli, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation

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

    ax.set_yticks(numpy.arange(len(data.keys()))+1)
    ax.set_yticklabels(data.keys())

    gtm = '#FF8A87' # red
    ltm = '#92A7FF' # blue

    xlim = ax.get_xlim()[1]
    if xlims:
        ax.set_xlim(xlims)
        xlim = xlims[1]

    draw_qs = True

    for i, k, p in zip(range(0, len(data)), data, r['boxes']):
        p.set_facecolor(col)

    if title:
        ax.set_title(title, fontsize=6)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)
    plot.close(fig)

'''
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
'''
