#!/bin/bash
ulimit -s unlimited
start_time=`date +%s`
count=1
result=log/`date +%Y%m%d_%H-%M-%S`.log
echo "---GRSA---"
for lambda in 5
do
    for T in  2 3 4 5 6 7 8
    do

    # range_size <= label_size
    job_start=`date +%s`
    ./grsa input/tsukuba_ output/tsukuba_${lambda}_${T}.bmp ${T} 16 $lambda >>  ${result}
    job_end=`date +%s`
    time=$((job_end - job_start));
    count=`expr $count + 1`
    echo "tsukuba, lambda=${lambda} [${time}s]";
    done
done


# # for T in 2 3 4 5 6 7 8
# # do
# #     for lambda in 1 2 3 4 5 6 7
# #     do
#
# #         # range_size <= label_size
# #         job_start=`date +%s`
# #         ./grsa input/venus_ output/venus_${lambda}_${T}.bmp $T 8 $lambda >>  ${result}
# #         job_end=`date +%s`
# #         time=$((job_end - job_start));
# #         count=`expr $count + 1`
# #         echo "venus T=${T}, lambda=${lambda} [${time}s]";
# #     done
# # done
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
time=$((end_time -start_time));
rm temp.txt

echo "@trsk_1st 全ての処理が完了しましたっ(GRSA)! 総所要時間[${time}s]" | tw --pipe --user="trsk_1st"

git add ${result}
git commit -m "job_${result}"
git push origin master

echo "------"
