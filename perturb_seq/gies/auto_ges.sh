#!/bin/sh

for file in `ls ./data/*.rda`
do
	for k in 0 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150
	do
		file=`basename $file`
		cat ges.R | R --slave --args "./data/"${file} ${k} ${k} ${k}
	done
done
