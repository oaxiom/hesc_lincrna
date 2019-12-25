
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
        destination = 'novel_transcript_%s_%s_%s' % (part3[sum[1]], part2[sum[0]], part4[sum[3]], part5[sum[2]])
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
    cdsl = gencode['cds_loc']['left']
    cdsr = gencode['cds_loc']['right']

    global_left = gencode['loc']['left']
    global_right = gencode['loc']['right']

    splice_sites = []
    newcdsl = 0 ; newcdsr = 0

    tlen = get_transcript_length(gencode)

    # convert all of the positions into the spliced locations, on the + strand;
    if gencode['strand'] == '+':
        left = 0 # by definition;
        rite = tlen


    elif gencode['strand'] == '-':
        left = 0 # by definition;
        rite = tlen

    '''
    exonStarts = gencode['exonStarts']
    exonEnds = gencode['exonEnds']

    for exon in zip(gencode['exonStarts'], gencode['exonEnds']):
        print(exon)
        if cdsl >= exon[0] and cdsl <= exon[1]:
            newcdsl = tlength + (cdsl-exon[0])

        if cdsr >= exon[0] and cdsr <= exon[1]:
            newcdsr = tlength + (cdsr-exon[0])

        tlength += (exon[1]-exon[0])
        currpos = exon[1]
        splice_sites.append(tlength)
    '''
        #print(tlength, cdsl, cdsr, splice, newcdsl, newcdsr, cdsl >= splice[0] and cdsl <= splice[1], cdsr >= splice[0] and cdsr <= splice[1])

    #if gencode and cdsl != cdsr: # if cdsl= cdsr then it is non-coding;
    cdsl = newcdsl
    cdsr = newcdsr

    splice_sites = splice_sites[:-1] # last one is the termination;

    #ts = gene['loc']['left']
    #te = gene['loc']['right']

    return 0, tlength, cdsl, cdsr, splice_sites

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
    currpos = gencode['loc']['left'] # transcript position in genomic coords
    splice_sites = []
    newcdsl = 0 ; newcdsr = 0
    for splice in zip(gencode['exonStarts'], gencode['exonEnds']):
        tlength += (splice[1]-splice[0])
        currpos = splice[1]
        splice_sites.append(tlength)

    splice_sites = splice_sites[:-1] # last one is the termination;

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
