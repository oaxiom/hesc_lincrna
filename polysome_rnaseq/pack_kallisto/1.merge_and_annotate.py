#! /usr/bin/env python3

import sys, os, glob, gzip, math
from glbase3 import *

all = None

conds = []

def pretify_sample_names(cond_names):
    newc = []
    for c in cond_names:
        t = c.replace('hesc_', 'hESC ').replace('_rp1', '').replace(' rp2', '#2')
        t = t.replace('nuc', 'Nuclear').replace('cyto', 'Cytoplasm').replace('monosome', 'Monosome')
        t = t.replace('poly_high', 'Polysome High').replace('poly_low', 'Polysome Low')
        t = t.replace('_', ' ')
        newc.append(t)
    return newc

# Failed QC:
bad_samples = [
    ]

all = None

conds = []
all_files = glob.glob("../kallisto/*/abundance.tsv")
all_files.sort()

for f in all_files:
    head = os.path.split(f)[0].split('/')[2]

    if ".rp" in head:
        head = "_".join(head.split(".")[0:2])
    else:
        head = head.split(".")[0]

    if head in bad_samples:
        continue

    tbase = head
    n = 1
    while tbase in conds: # Make sure all condition names are unique.
        tbase = "%s_%s" % (tbase, n)
    base = tbase

    conds.append(base)

    print("...", base)

    rsem = []
    oh = open(f, "r")
    for line in oh:
        ll = line.strip().split("\t")
        if not "target_id" in ll[0]:
            rsem.append({"transcript_id": ll[0], base: float(ll[4])})
    oh.close()

    if not all:
        all = rsem
    else:
        for idx, g in enumerate(all):
            if g["transcript_id"] == rsem[idx]["transcript_id"]:
                g[base] = rsem[idx][base]
            else:
                print("Oh dear, died on sample: %s" % f)
                sys.exit()

cond_names = conds

custom_ann = glload('../../transcript_assembly/packed/all_genes.glb')
expn = expression(loadable_list=all, expn=cond_names)

print(expn)
print(custom_ann)

expn = custom_ann.map(genelist=expn, key='transcript_id')
expn = expn.getColumns(['ensg', 'enst', 'name', 'transcript_id', 'transcript_class', 'coding', 'expression'])

expn.saveTSV("kall_tpm-unmerged.tsv")
expn.save("kall_tpm-unmerged.glb") # Only need for QC purposes

print(expn.getConditionNames())

# Merge:
expn = expn.mean_replicates(
    ['hesc_cyto_rp1', 'hesc_cyto_rp2'],
    ['hesc_nuc_rp1', 'hesc_nuc_rp2'],
    ['hesc_monosome_rp1', 'hesc_monosome_rp2'],
    ['hesc_poly_low_rp1', 'hesc_poly_low_rp2',],
    ['hesc_poly_high_rp1', 'hesc_poly_high_rp2'],
    output_pears="Pearson_correlations.tsv",
    pearson_hist="Pearson_hist.png",
    threshold=0.8)

expn.strip_errs()

expn.setConditionNames(pretify_sample_names(expn.getConditionNames()))

print(expn.getConditionNames())

expn = expn.sliceConditions(['hESC Nuclear', 'hESC Cytoplasm', 'hESC Monosome', 'hESC Polysome Low', 'hESC Polysome High', ])

expn.saveTSV("kall_tpm.tsv")
expn.save("kall_tpm.glb")

# convert the TPMs to normalised ratios:

expn_poly = expn.sliceConditions(['hESC Cytoplasm', 'hESC Monosome', 'hESC Polysome Low', 'hESC Polysome High', ])
newe = []
for e in expn_poly:
    s = e['conditions'][0]
    if s > 0.01: # one or two nan. This seq is not as deep as our full dataset, so some transcripts become unreliable.
        e['conditions'] = [i/s for i in e['conditions']]
        newe.append(e)
expn_poly.load_list(newe, cond_names=expn_poly.getConditionNames())
expn_poly.saveTSV("kall_poly_ratios.tsv")
expn_poly.save("kall_poly_ratios.glb")

expn_nuccyt = expn.sliceConditions(['hESC Nuclear', 'hESC Cytoplasm'])
newe = []
for e in expn_nuccyt:
    s = sum(e['conditions'])
    if e['conditions'][0] > 0.01: # one or two nan. This seq is not as deep as our full dataset, so some transcripts become unreliable.
        e['conditions'] = [math.log2(e['conditions'][1]/e['conditions'][0]), ]
        # i.e. lncATALAS RIC
        newe.append(e)
expn_nuccyt.load_list(newe, cond_names=['RIC'])
expn_nuccyt.saveTSV("kall_nuc_cyt_ratios_RIC.tsv")
expn_nuccyt.save("kall_nuc_cyt_ratios_RIC.glb")
