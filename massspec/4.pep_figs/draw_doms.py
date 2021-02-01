
from glbase3 import glload, genelist
import draw_domains_share

res = glload('../results_gene.glb')

doms = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
CDSs = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
CDSs = {i['transcript_id']: i for i in CDSs}
dfam = genelist('../../te_discovery/dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

draw = 'pdf'

for n, gene in enumerate(res):
    _class = gene['class']
    gene = doms.get(key='transcript_id', value=gene['transcript_id']) # By definition

    if not gene:
        continue

    gene = gene[0]

    peptides = res.get(value=gene['transcript_id'], key='transcript_id', mode='greedy')
    if not peptides:
        continue

    print(gene['transcript_id'])

    if gene['transcript_id'] in CDSs:
        gene['cds_local_locs'] = CDSs[gene['transcript_id']]['cds_local_locs']

    peptides = [(i['mrna_left'], i['mrna_right'], i['peptide_string'], i['insideTE']) for i in peptides]

    draw_domains_share.draw_domain(gene,
        '%s/%s.%s.%s.%s.%s' % (draw, _class, gene['name'], gene['transcript_id'], gene['enst'], draw),
        dfam,
        peptides
        )


