
gunzip -c ../gtf/current_gtf.gtf.gz | sed 's#^#chr#g' >current_gtf.fixed_contigs.gtf
