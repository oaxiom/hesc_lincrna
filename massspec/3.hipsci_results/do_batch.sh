
for f in  ../centroid_mzml/*.cwt.mzML
do
    bf=`basename $f | sed -r 's#.cwt.mzML##g' `
    tt=`echo $bf.tsv`
    
    if [ ! -f $tt ] && [ ! -f $tt.gz ]
    then
        echo ... $tt
        qsub -N ms.${bf} -v in=${f},out=${bf} msgf_search.sh
        sleep 1
    fi
done

