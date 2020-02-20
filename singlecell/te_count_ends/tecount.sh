#PBS -l nodes=1:ppn=2,mem=32gb
#PBS -j oe
#PBS -o ${out}.txt
#PBS -q batch
#PBS -V 
cd $PBS_O_WORKDIR

te_count -m custom -g 3ends_only.idx -i $inp -o ${out}.tsv --sc --se --strand -w ${white}
gzip ${out}.tsv

