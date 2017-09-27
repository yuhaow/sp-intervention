#!/usr/bin/env bash

for file in `ls ./greedy_sp_data_noisy/*.rda`
do
	for method in hsic.gamma
	do
		for alpha in 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6
		do
			file=`basename $file`
			nohup cat sp_restart.R | R --slave --args "./greedy_sp_data_noisy/"${file} ${alpha} ${method} &

		done
	done
done

#   hsic.clust
