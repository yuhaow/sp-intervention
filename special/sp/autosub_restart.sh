#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	for alpha in 0.0001 0.0005 0.001 0.005 0.01 0.05
	do
		file=`basename $file`
		echo "\$R --no-save --args ../data/"${file} ${alpha}" < sp_restart.R" > tmp.sh
		cat template.sh tmp.sh > scripts/${file}.${alpha}.sh
		rm tmp.sh
		sbatch scripts/${file}.${alpha}.sh
	done
done
