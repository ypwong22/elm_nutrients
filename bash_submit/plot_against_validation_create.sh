create=1
if [ ${create} -eq 1 ]
then
    for i in {0..9}
    do
	for j in {0..9}
	do
	    cat ../plot_against_validation.py | sed s/REPLACE1/${i}/ | sed s/REPLACE2/${j}/ > ../plot_against_validation_${i}_${j}.py
	    cat plot_against_validation_template.sh | sed s/REPLACE1/${i}/ | sed s/REPLACE2/${j}/ > plot_against_validation_${i}_${j}.sh
	done
    done
fi


clean=0
if [ ${clean} -eq 1 ]
then
    for i in {0..9}
    do
	for j in {0..9}
	do
	    rm plot_against_validation_${i}_${j}.sh
	done
    done
fi


submit=1
if [ ${submit} -eq 1 ]
then
    for i in {0..9}
    do
	for j in {0..9}
	do
	    qsub plot_against_validation_${i}_${j}.sh
	done
    done
fi
