
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
