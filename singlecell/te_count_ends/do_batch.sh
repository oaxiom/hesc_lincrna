
for f in  ../starsolo.v1/*.Aligned.sortedByCoord.out.bam
do
    root=`basename $f`
    bf=`echo $root | sed -r 's/.Aligned.sortedByCoord.out.bam//g' | sed -r 's/.gz//g'`
    p1=`echo $root`
    tt=`echo $bf.tsv.gz`
    if [ ! -f $tt ]
    then
        echo PE ... $tt 
        qsub -N tec.$bf -v inp=$f,out=$bf,white=version1.txt tecount.sh
        sleep 1
    fi
done

for f in  ../starsolo.v2/*.Aligned.sortedByCoord.out.bam
do
    root=`basename $f`
    bf=`echo $root | sed -r 's/.Aligned.sortedByCoord.out.bam//g' | sed -r 's/.gz//g'`
    p1=`echo $root`
    tt=`echo $bf.tsv.gz`
    if [ ! -f $tt ]
    then
        echo PE ... $tt 
        qsub -N tec.$bf -v inp=$f,out=$bf,white=version2.txt tecount.sh
        sleep 1
    fi
done

