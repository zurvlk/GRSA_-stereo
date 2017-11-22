#!/bin/bash

start_time=`date +%s`
count=1
result=log/`date +%Y%m%d_%H-%M-%S`.log
input_file="input/barbara_noise32.bmp"
for func in 1
do
    for label_size in 16 32 64
    do
        for range_size in  4 8 16 32 64
        do
            if [ $range_size -le $label_size ]
            then
                # range_size <= label_size
                job_start=`date +%s`
                ./grsa ${input_file} output/output_${func}_${label_size}_${range_size}.bmp $label_size $range_size $func | tee temp.txt >>  ${result}
                job_end=`date +%s`
                time=$((job_end - job_start));
                echo "@zurvlk 処理${count}が完了しましたっ!" > notification.txt

                if [ $func -eq 0 ]
                then 
                    echo "Vpq(fp,fq) = |fp - fq|" >> notification.txt
                else
                    echo "Vpq(fp,fq) = (fp - fq)^2" >> notification.txt
                fi

                # echo "input_file: ${input_file}"
                echo "label_size: ${label_size} range_size: ${range_size}" >> notification.txt
                cat temp.txt | grep before >> notification.txt
                cat temp.txt | grep after >> notification.txt

                cat notification.txt | tw --pipe --user="trsk_1st"
                count=`expr $count + 1`
            fi
        done
    done
done
end_time=`date +%s`
time=$((end_time -start_time));
rm temp.txt

echo "@zurvlk 全ての処理が完了しましたっ! 総所要時間[${time}s]" | tw --pipe --user="trsk_1st"

git add ${result}
git commit -m "job_${result}"
git push origin master
