#!/bin/bash
ulimit -s unlimited
start_time=`date +%s`
count=1
result=`date +%Y%m%d_%H-%M-%S`.log
echo "---GRSA---"
for dterm in 0 1
do
    for lambda in 1 2 3 4 5 
    do
        for T in  2 3 4 5 6 7 8
        do
        # range_size <= label_size
            job_start=`date +%s`
            ./grsa input/tsukuba_ output/tsukuba_${dterm}_${lambda}_${T}.bmp ${T} 16 $lambda ${dterm}>>  log/${dterm}_lambda_${lambda}_tsukuba${result}
            job_end=`date +%s`
            time=$((job_end - job_start));
            count=`expr $count + 1`
            echo "tsukuba, lambda=${lambda} [${time}s]";
        done
        git add log/${dterm}_lambda_${lambda}_tsukuba${result}
        git commit -m "job_${result}"
        git push origin master
    done
done

for dterm in 0 1
do
    for lambda in 1 2 3 4 5
    do
        for T in 2 3 4 5 6 7 8
        do
        # range_size <= label_size
            job_start=`date +%s`
            ./grsa input/venus_ output/venus_${dterm}_${lambda}_${T}.bmp $T 8 $lambda ${dterm} >> log/${dterm}_lambda_${lambda}_venus${result}
            job_end=`date +%s`
            time=$((job_end - job_start));
            count=`expr $count + 1`
            echo "venus T=${T}, lambda=${lambda} [${time}s]";
        done
        git add log/${dterm}_lambda_${lambda}_venus${result}
        git commit -m "job_${result}"
        git push origin master
    done
done

for dterm in 0 1
do
    for lambda in 1 2 3 4 5
    do
        for T in 2 3 4 5 6 7 8
        do
        # range_size <= label_size
            job_start=`date +%s`
            ./grsa input/teddy_ output/teddy__${dterm}_${lambda}_${T}.bmp $T 4 $lambda ${dterm} >> log/${dterm}_lambda_${lambda}_teddy${result}
            job_end=`date +%s`
            time=$((job_end - job_start));
            count=`expr $count + 1`
            echo "teddy T=${T}, lambda=${lambda} [${time}s]";
        done
        git add log/${dterm}_lambda_${lambda}_teddy${result}
        git commit -m "job_${result}"
        git push origin master
    done
done


#
# for T in  5 6 7 8
# do
#     for lambda in 1 2 3 4 5 6 7
#     do
#         job_start=`date +%s`
#         ./grsa input/teddy_ output/teddy_${lambda}_${T}.bmp $T 4 $lambda >>  ${result}
#         job_end=`date +%s`
#         time=$((job_end - job_start));
#         count=`expr $count + 1`
#         echo "teddy T=${T}, lambda=${lambda} [${time}s]";
#     done
# done

end_time=`date +%s`
time=$((end_time - start_time));
rm temp.txt

echo "@trsk_1st 全ての処理が完了しましたっ(GRSA)! 総所要時間[${time}s]" | tw --pipe --user="trsk_1st"
echo "------"
