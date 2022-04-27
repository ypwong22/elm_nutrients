create=0
if [ ${create} -eq 1 ]
then
    for i in 04 06 07 08 10 11 13 16 17 19 20
    do
	cat plot_ts_multiruns_template.sh | sed s/REPLACE/${i}/ > plot_ts_multiruns_${i}.sh
    done
fi


clean=1
if [ ${clean} -eq 1 ]
then
    for i in 04 06 07 08 10 11 13 16 17 19 20
    do
	rm plot_ts_multiruns_${i}.sh
    done
fi


submit=0
if [ ${submit} -eq 1 ]
then
    for i in 04 06 07 08 10 11 13 16 17 19 20
    do
	qsub plot_ts_multiruns_${i}.sh
    done
fi
