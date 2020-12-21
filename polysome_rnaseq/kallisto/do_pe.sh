#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o ${out}.results.txt
#PBS -V 
cd $PBS_O_WORKDIR

#bamopts='--genomebam -g /share/apps/genomics/genome/mm10/ensembl/mm10.ensemblv92.nopsuedo.gtf -c /share/apps/genomics/genome/mm10/mm10.chromSizes' 
index='/data3/andrew/telncrna/gtf/kallisto_index/current_gtf.kall.idx'

kallisto quant -t 1 -b 10 --bias --plaintext -i $index -o ${out} ${p1} ${p2}  


