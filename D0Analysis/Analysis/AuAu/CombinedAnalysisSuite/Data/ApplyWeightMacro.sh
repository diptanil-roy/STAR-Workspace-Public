#!/bin/bash

# All D0 pT [1-10]

[[ -z "$1" ]] && { echo "Enter The Name For the sWeights Directory As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Name For the Histograms Directory As Argument #2" ; exit 1; }

cp applyweights.C applyweights_runner.C
sed -i'.bak' 's/applyweights/applyweights_runner/g' applyweights_runner.C
sed -i'.bak' "s/sWeights/$1/g" applyweights_runner.C

root -l -b -q applyweights_runner.C


for pt in {1..5}
do

	# # All D0 pT [1-pt]

	# cp applyweights.C applyweights_runner.C
	# sed -i'.bak' 's/applyweights/applyweights_runner/g' applyweights_runner.C
	# sed -i'.bak' "s/sWeights/$1/g" applyweights_runner.C
	# sed -i'.bak' "s/double pthigh = 10/double pthigh = $pt/" applyweights_runner.C

	# root -l -b -q applyweights_runner.C

	# All D0 pT [pt-10]

	cp applyweights.C applyweights_runner.C
	sed -i'.bak' 's/applyweights/applyweights_runner/g' applyweights_runner.C
	sed -i'.bak' "s/sWeights/$1/g" applyweights_runner.C
	sed -i'.bak' "s/double ptlow = 1/double ptlow = $pt/" applyweights_runner.C

	root -l -b -q applyweights_runner.C
	root -l -b -q 'applyweights_runner.C(1)'
	root -l -b -q 'applyweights_runner.C(-1)'
	root -l -b -q 'applyweights_runner.C(2)'
	root -l -b -q 'applyweights_runner.C(-2)'

done

mkdir -p $2
mv Histograms* $2/.

rm applyweights_runner.C
rm applyweights_runner.C.bak