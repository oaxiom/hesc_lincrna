

for f in *.fa
do
    o=`basename $f | sed 's#.fa##g'` 

    qsub -v inp=$f,out=$o hmmsearch.pbs
done
