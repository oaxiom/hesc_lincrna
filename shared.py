
col_keys = {
    'DNA:.': 'royalblue',
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

    'LTR:.': 'coral',
    'LTR:ERV1': 'tomato',
    'LTR:ERVK': 'orangered',
    'LTR:ERVL': 'darkred',
    'LTR:ERVL-MaLR': 'brown',
    'LTR:Gypsy': 'indianred',
    'LINE:RTE-X': 'tomato',

    'RC:Helitron': 'greenyellow',

    'Retroposon:SVA': 'orange', # To get new cols;

    'Satellite:.': 'olive',
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

def get_cols(labels):
    cols = []
    for e in labels:
        cols.append(col_keys[e])
    return(cols)

part1 = {'~)': 'novel_isoform',
    '=)': 'known_isoform'}

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
    if 'HSC' in name:
        destination = 'unknown_transcript_%s_%s_%s' % (part3[sum[1]], part2[sum[0]], part4[sum[2]])
        alpha = '.'

    else:
        destination = '%s_%s_%s_%s' % (part1[sum[3]], part3[sum[1]], part2[sum[0]], part4[sum[2]])

        alpha = name[0]

    return(destination, alpha)

if __name__ == '__main__':
    import matplotlib.pyplot as plot

    # legend:
    fig, ax=plot.subplots(figsize=(4,6))
    ax.set_xlim(0, 30)
    ax.set_axis_off()

    for i, name in enumerate(reversed(list(col_keys.keys()))):
        y = i * 10

        ax.text(4, y, name, fontsize=6,
                horizontalalignment='left',
                verticalalignment='center')

        ax.hlines(y, 1, 2, color=col_keys[name], linewidth=5)
    fig.savefig('legend.png')
    fig.savefig('legend.svg')

    print('\n'.join(list(col_keys.keys())))
