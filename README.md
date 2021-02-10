# hesc_lincrna

Transposable elements (TEs) occupy nearly 50% of mammalian genomes and are both potential dangers to genome stability and functional genetic elements. TEs can be expressed and exonised as part of a transcript, however, their full contribution to the transcript splicing remains unresolved. Here, guided by long and short read sequencing of RNAs, we show that 26% of coding and 65% of noncoding transcripts of human pluripotent stem cells (hPSCs) contain TEs. Different TE families have unique integration patterns with diverse consequences to transcript expression and function. These effects are at least partly mediated by TE-type specific changes in transcript cellular localization, alterations in RNA stability, and TE-specific binding of RBPs. We identify hPSC-specific splicing of endogenous retroviruses (ERVs) as well as LINE L1 elements into protein coding genes that generate TE-derived peptides. Finally, single cell RNA-seq reveals that proliferating hPSCs are dominated by ERV-containing transcripts, and subpopulations express SINE or LINE-containing transcripts. Overall, we demonstrate that TEs are widely found inside the transcriptome of pluripotent cells, are modulated by several cellular mechanisms, and lead to novel transcripts and peptides.   

Preprint:

https://www.biorxiv.org/content/10.1101/2020.07.26.220608v1

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