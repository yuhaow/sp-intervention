#!/usr/bin/env bash

for file in `ls ../data/*.rda`
do
	file=`basename $file`
	echo "\$R --no-save --args ../data/"${file}" < ges.R" > tmp.sh
	cat template.sh tmp.sh > scripts/${file}.sh
	rm tmp.sh
	sbatch scripts/${file}.sh
done
