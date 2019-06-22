
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

part1 = {'~)': 'variant_isoform',
    '=)': 'known_isoform',
    '!)': 'unknown_isoform',}

part2 = {'(ME': 'multi_exon',
    '(SE': 'single_exon'}

part3 = {'C': 'coding',
    'NC': 'noncoding',
    'U': 'unknown'}

part4 = {'SR': 'short_read',
    'SR+LR': 'long_read', # Takes preference for purposes of prrof
    'LR': 'long_read'}

def classify_transcript(name):
    sum = name.split(' ')[1].split(';')
    if 'HSC' in name: # and also ';!)'
        destination = 'novel_transcript_%s_%s_%s' % (part3[sum[1]], part2[sum[0]], part4[sum[2]])
        alpha = '.'

    else:
        destination = '%s_%s_%s_%s' % (part1[sum[3]], part3[sum[1]], part2[sum[0]], part4[sum[2]])

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
    # This is wrong if the transcrip is ~
    cdsl = gencode['cds_loc']['left']
    cdsr = gencode['cds_loc']['right']
    #print(gencode['enst'], gencode['transcript_id'], gencode['name'], gencode['strand'], gencode['loc'], gencode['cds_loc'], gencode['exonStarts'], gencode['exonEnds'])

    # Work out the mRNA length from the gene length and the exons and Es:
    tlength = 0
    currpos = gencode['loc']['left'] # transcript position in genomic coords
    splice_sites = []
    newcdsl = 0 ; newcdsr = 0
    for splice in zip(gencode['exonStarts'], gencode['exonEnds']):
        if cdsl >= splice[0] and cdsl <= splice[1]:
            newcdsl = tlength + (cdsl-splice[0])
        if cdsr >= splice[0] and cdsr <= splice[1]:
            newcdsr = tlength + (cdsr-splice[0])

        tlength += (splice[1]-splice[0])
        currpos = splice[1]
        splice_sites.append(tlength)

        #print(tlength, cdsl, cdsr, splice, newcdsl, newcdsr, cdsl >= splice[0] and cdsl <= splice[1], cdsr >= splice[0] and cdsr <= splice[1])

    #if gencode and cdsl != cdsr: # if cdsl= cdsr then it is non-coding;
    cdsl = newcdsl
    cdsr = newcdsr

    splice_sites = splice_sites[:-1] # last one is the termination;

    #ts = gene['loc']['left']
    #te = gene['loc']['right']

    return 0, tlength, cdsl, cdsr, splice_sites

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
