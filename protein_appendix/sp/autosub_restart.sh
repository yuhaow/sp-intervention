#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for alpha in 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${alpha}" < sp_restart.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.sh
		rm tmp.sh
		sbatch scripts/${file}.sh
	done
done
