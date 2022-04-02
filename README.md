# hesc_lincrna

Transposable elements (TEs) occupy nearly 40% of mammalian genomes and, whilst most are fragmentary and no longer capable of transposition, they can nevertheless contribute to cell function. TEs within genes transcribed by RNA polymerase II can be copied as parts of primary transcripts; however, their full contribution to mature transcript sequences remains unresolved. Here, using long and short read (LR and SR) RNA sequencing data, we show that 26% of coding and 65% of noncoding transcripts in human pluripotent stem cells (hPSCs) contain TE-derived sequences. Different TE families are incorporated into RNAs in unique patterns, with consequences to transcript structure and function. The presence of TE sequences within a transcript is correlated with TE-type specific changes in its subcellular distribution, alterations in steady-state levels and half-life, and differential association with RNA Binding Proteins (RBPs). We identify hPSC-specific incorporation of endogenous retroviruses (ERVs) and LINE:L1 into protein-coding mRNAs, which generate TE sequence-derived peptides. Finally, single cell RNA-seq reveals that hPSCs express ERV-containing transcripts, whilst differentiating subpopulations lack ERVs and express SINE and LINE-containing transcripts. Overall, our comprehensive analysis demonstrates that the incorporation of TE sequences into the RNAs of hPSCs is more widespread and has a greater impact than previously appreciated.

Published online in Nucleic Acids Research:

https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkab710

Required Installation:

1. Python >=3.6
2. Numpy, scipy, sklearn, matplotlib
3. scanpy, 
4. glbase3
5. RSEM and hg38 genome
6. chipFish (Optional)

Code run order:

1. Execute gencode/download_gencode.sh
2. Execute genome_repeats/get_rmsk.sh
3. Execute transcript_assembly/get_CDS/download_data.sh
5. Pack and annotate the GTF: transcript_assembly/packed/pack_gtf.py
6. Fix the GTF contig names: transcript_assembly/fasta/fix_gtf_contig_names.sh
7. Get the FASTA using RSEM (qsub, or just execute): transcript_assembly/fasta/gtf_to_fasta.pbs
8. Build the TE indeces: te_discovery/te_transcripts/get_te_containing_transcripts.py

From there, it should be possible to pretty much run the py scripts in any order. 
Some may give errors, if they require certain other processed files to exist, but the above is the only critical order.
