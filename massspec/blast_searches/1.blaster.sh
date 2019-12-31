
blastp -ungapped -comp_based_stats F -query ../get_pep_fastas/fasta/table_new_ATG.fasta -evalue 1e-20 -outfmt 6 -num_threads 2 -db ../blast_db/gencode.v32.pc_translations.fa -out blaster/table_new_ATG.tsv
blastp -ungapped -comp_based_stats F -query ../get_pep_fastas/fasta/table_new_ATG.fasta -evalue 1e-20 -outfmt 1 -num_threads 2 -db ../blast_db/gencode.v32.pc_translations.fa -out blaster/table_new_ATG-aligns.txt


