#!/bin/sh

for file in `ls ./greedy_sp_data/*.rda`
do

	file=`basename $file`
	cat sp_gauss.R | R --slave --args "./greedy_sp_data/"${file}

done
