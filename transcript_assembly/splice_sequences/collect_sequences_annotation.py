
import sys, os, gzip
from glbase3 import genome, utils

hg38 = genome()
hg38.bindSequence(os.path.expanduser('~/hg38/seq/'))

oh = gzip.open('../gtf/current_gtf.gtf.gz', 'rt') # you need to download this one

p5 = {'matching':
        {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0},
    'variant': {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0},
    'novel': {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0}
    }


p3 = {'matching':
        {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0},
    'variant': {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0},
    'novel': {'AA': 0, 'AC': 0, 'AG': 0, 'AT': 0,
        'CA': 0, 'CC': 0, 'CG': 0, 'CT': 0,
        'GA': 0, 'GC': 0, 'GG': 0, 'GT': 0,
        'TA': 0, 'TC': 0, 'TT': 0, 'TG': 0,
        'NN': 0}
    }

tss = []
tts = [] # not accurately named

done = {'matching': 0, 'variant': 0, 'novel': 0}
skipped = 0
trans = None
transcript_types_seen = set([])
transcript_types_dont_keep = set([
    'miRNA', 'snRNA', 'snoRNA', 'Mt_tRNA', 'Mt_rRNA', 'TEC', 'scRNA', 'scaRNA','vaultRNA',
    'ribozyme', 'transcribed_unprocessed_pseudogene', 'unprocessed_pseudogene',
    'IG_pseudogene', 'IG_J_gene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'TR_D_gene',
    'IG_V_gene', 'IG_C_pseudogene',  'TR_J_gene','TR_V_pseudogene','sRNA', 'TR_V_gene', 'IG_D_gene', 'IG_C_gene', 'IG_V_pseudogene', 'TR_C_gene',])

for idx, line in enumerate(oh):
    if '#' in line[0]:
        continue
    line = line.strip().split('\t')
    #print(line)

    if line[2] == 'start_codon': continue
    if line[2] == 'stop_codon': continue
    if line[2] == 'CDS': continue
    if line[2] == 'UTR': continue
    if line[2] == 'gene': continue

    if line[2] == 'transcript':
        if line[6] == '+':
            tss = int(line[3])
            tts = int(line[4])
        elif line[6] == '-':
            tss = int(line[4])
            tts = int(line[3])
        continue

    # collect exon seqs;
    gtf_dec = {}
    for item in line[8].split(';'):
        if item:
            item = item.strip(' ').replace('"', '').split(' ')
            gtf_dec[item[0]] = item[1]

    if 'gene_type' in gtf_dec and gtf_dec['gene_type'] in transcript_types_dont_keep:
        continue
    if 'transcript_type' in gtf_dec and gtf_dec['transcript_type'] in transcript_types_dont_keep:
        continue

    l = int(line[3])
    r = int(line[4])
    c = 'chr{}'.format(line[0])
    k = gtf_dec['transcript_class']

    if line[6] == '+':
        # check the 5' of the 'exon' is not the TSS:
        if l != tss:
            seq = hg38.getSequence('{}:{}-{}'.format(c, l-2, l-1)).upper() # ....NN|RNA
            p3[k][seq] += 1

        if r != tts:
            seq = hg38.getSequence('{}:{}-{}'.format(c, r+1, r+2)).upper() # RNA|NN...
            p5[k][seq] += 1

    elif line[6] == '-':
        # check the 5' of the 'exon' is not the TSS:
        if l != tss:
            seq = utils.rc(hg38.getSequence('{}:{}-{}'.format(c, l-2, l-1)).upper()) # ....NN|RNA
            p5[k][seq] += 1

        if r != tts:
            seq = utils.rc(hg38.getSequence('{}:{}-{}'.format(c, r+1, r+2)).upper()) # RNA|NN...
            p3[k][seq] += 1

    done[k] += 1

    if idx % 1e4 == 0:
        print('{:,}'.format(idx))
        #break

print(p5)
print(p3)

for k in p5:
    oh = open('hpsc_transcriptome_{}_results.tsv'.format(k), 'wt')
    oh.write("sequence\t5'\t5'%\t3'\t3'%\n")
    for seq in sorted(p5[k].keys()):
        oh.write('{}\t{}\t{:.2f}%\t{}\t{:.2f}%\n'.format(seq, p5[k][seq], p5[k][seq]/done[k]*100, p3[k][seq], p3[k][seq]/done[k]*100))
    oh.close()
