#! /usr/bin/env python3

import sys, os, glob, gzip
from glbase3 import *

sys.path.append('../')
import shared

all = None

conds = []

# Failed QC:
bad_samples = [
    #'Hs_esc_naive_win1rj_rp2',
    #'Hs_esc_naive_win1rj_rp1',
    #'Hs_esc_naive_wibr3rjtgcl16_rp1',
    #'Hs_esc_naive_wibr3rjtgcl12_rp1',
    #'Hs_esc_naive_wibr2rj_rp1',
    ]

all = None

conds = []
all_files = glob.glob("data/*/*.tsv")
all_files.sort()

for f in all_files:
    head = os.path.split(f)[0].split('/')[1]

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

expn = expn.filter_low_expressed(5, 2) # works out to about ~30 mapped tags
expn = custom_ann.map(genelist=expn, key='transcript_id')
#expn.coerce(int) #?!?!
expn.saveTSV("kall_tpm-unmerged.tsv")
expn.save("kall_tpm-unmerged.glb") # Only need for QC purposes

print(expn.getConditionNames())

# Merge:
expn = expn.mean_replicates(
    ['Hs_esc_primed_H1_day0_rp1', 'Hs_esc_primed_H1_day0_rp2'],
    ['Hs_tsc_rp1', 'Hs_tsc_rp2', 'Hs_tsc_rp3', 'Hs_tsc_rp4', 'Hs_tsc_rp5', 'Hs_tsc_rp6'],
    ['Hs_hesc_hde_d0_rp1', 'Hs_hesc_hde_d0_rp2', 'Hs_hesc_hde_d0_rp3'],
    ['Hs_hesc_hde_d5_rp1', 'Hs_hesc_hde_d5_rp2', 'Hs_hesc_hde_d5_rp3'],
    output_pears="Pearson_correlations.tsv",
    pearson_hist="Pearson_hist.png",
    threshold=0.8)

expn = expn.filter_low_expressed(5, 1) # works out to about ~30 mapped tags
expn.setConditionNames(shared.pretify_sample_names(expn.getConditionNames()))
expn.saveTSV("kall_tpm.tsv")
expn.save("kall_tpm.glb") # Only need for QC purposes


