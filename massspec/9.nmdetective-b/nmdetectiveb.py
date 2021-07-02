
import sys, os
from glbase3 import *
sys.path.append('../../')
import shared

def NMDetective_B_score(inlastexon, orflength, exonlength, within_50nt_of_lastEJ):
    if inlastexon:
        return 0.0, 'Last exon'

    if orflength < 150: # Distance from STOP to START
        return 0.12, 'Start-proximal'

    if exonlength >407: # exon length the STOP codon is in
        return 0.41, 'Long exon'

    if within_50nt_of_lastEJ: # EJC will block NMD in STOP within 50 nt.
        return 0.2, '50 nt rule'

    return 0.65, 'Trigger NMD'

all_transcripts = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
all_exons = glload('../../transcript_assembly/packed/all_genes.glb')

all_transcripts = all_transcripts.map(key='transcript_id', genelist=all_exons)
all_expressed = all_transcripts.getColumns(['enst', 'transcript_id'])
all_coding = glload('../../gencode/hg38_gencode_v32.pc.glb').map(genelist=all_transcripts, key='enst')

matching_coding = all_coding.get('transcript_class', 'matching')
variant_coding = all_coding.get('transcript_class', 'variant')

has_peptide_hit = glload('../results_gene.glb').removeDuplicates('transcript_id') # ones with a peptide hit.
has_peptide_hit = has_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

no_peptide_hit = glload('../2.blast_searches/super_table.glb').removeDuplicates('transcript_id')
no_peptide_hit = has_peptide_hit.map(genelist=no_peptide_hit, key='transcript_id', logic='notright')
no_peptide_hit = no_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

res = {
    'MS Hit': [],
    'No MS hit': [],
    'Variant': [],
    'Matching': [],
    }

data = {
    'Matching': matching_coding,
    'Variant': variant_coding,
    'MS Hit': has_peptide_hit,
    'No MS hit': no_peptide_hit,
    }

__notfound = 0
__nostop = 0

for k in data:
    for transcript in data[k]:
        strand = transcript['strand']

        cds_key_to_use = None
        if 'cds_loc' in transcript and transcript['cds_loc']:
            cds_key_to_use = 'cds_loc'
        elif 'cds_genome_loc' in transcript and transcript['cds_genome_loc']:
            cds_key_to_use = 'cds_genome_loc'
        else:
            cds_key_to_use = 'cds_local_to_genome'

        orflength = abs(transcript[cds_key_to_use]['right'] - transcript[cds_key_to_use]['left'])

        if transcript[cds_key_to_use]['right'] == 0:
            # Ignore misannotated GENCODE CDSs
            continue

        if strand == '+':
            STOP = transcript[cds_key_to_use]['right']
        else:
            STOP = transcript[cds_key_to_use]['left']

        inlastexon = False
        within_50nt_of_lastEJ = False

        if STOP == 0:
            __nostop += 1
            continue

        #print(strand, list(zip(transcript['exonStarts'], transcript['exonEnds'])))
        for exon_num, exons in enumerate(zip(transcript['exonStarts'], transcript['exonEnds'])):
            if STOP >= exons[0] and STOP <= exons[1]:
                if strand == '+':
                    d = exons[1] - STOP
                    if exon_num == len(transcript['exonEnds'])-1:
                        inlastexon = True
                        if d < 50:
                            within_50nt_of_lastEJ = True

                else:
                    d = STOP - exons[0]
                    if exon_num == 0:
                        inlastexon = True
                        if d < 50:
                            within_50nt_of_lastEJ = True

                exonlength = exons[1] - exons[0]

        nmd_score, nmd_class = NMDetective_B_score(inlastexon, orflength, exonlength, within_50nt_of_lastEJ)
        #print(nmd_score, transcript[cds_key_to_use], strand, inlastexon, orflength, exonlength, within_50nt_of_lastEJ)

        res[k].append(nmd_score)


shared.boxplots_simple('box_nmdb_estimate.pdf', res, None,
    #xlims=[0, 2500],
    col='#FF8A87',
    #vlines=[0, 55, 1000],
    showmeans=True
    )

shared.violin_simple('viol_nmdb_estimate.pdf', res, None,

    col='#FF8A87'
    )

shared.plots_simple('simp_nmdb_estimate.pdf', res, None,
    xlims=[-0.05, 0.55],
    vlines=[0, 0.25, 0.5],
    col='#FF8A87',
    )

