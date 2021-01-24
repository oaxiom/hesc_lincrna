#PBS -N msgfp.${out}
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -o ${out}.txt
#PBS -V 
cd $PBS_O_WORKDIR

opts1='-e 1' # Lys-C and Trypsin, Lys-C is a subset of Trypsin, so use Trypsin
opts2='-inst 1' # Was an orbitrap
#opts3='-t 0.5Da,2.5Da' # loose precision
opts3='-t 20ppm -ti -1,2' # high precision
opts4='-mod mods.mods'
java -Xmx3500M -jar ../../mgsfplus/MSGFPlus.jar -s ${in} -thread 30 -d all_masked_peptides.fa $opts1 $opts2 $opts3 $opts4 -ntt 2 -tda 1 -o ${out}.mzid 

# And cleanup
java -Xmx3500M -cp ../../mgsfplus/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ${out}.mzid -o ${out}.tsv
gzip ${out}.tsv
rm ${out}.mzid

#Example (high-precision): java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 1 -t 20ppm -ti -1,2 -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt
#Example (low-precision):  java -Xmx3500M -jar MSGFPlus.jar -s test.mzML -d IPI_human_3.79.fasta -inst 0 -t 0.5Da,2.5Da    -ntt 2 -tda 1 -o testMSGFPlus.mzid -mod Mods.txt

