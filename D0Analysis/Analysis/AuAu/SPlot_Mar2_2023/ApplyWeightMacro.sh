#!/bin/bash

# All D0 pT [1-10]

cp applyWeightsFor1GeV.C applyWeightsFor1GeV_runner.C
sed -i'.bak' 's/applyWeightsFor1GeV/applyWeightsFor1GeV_runner/g' applyWeightsFor1GeV_runner.C

root -l -b -q applyWeightsFor1GeV_runner.C


for pt in {2..5}
do

	# All D0 pT [1-pt]

	cp applyWeightsFor1GeV.C applyWeightsFor1GeV_runner.C
	sed -i'.bak' 's/applyWeightsFor1GeV/applyWeightsFor1GeV_runner/g' applyWeightsFor1GeV_runner.C
	sed -i'.bak' "s/double pthigh = 10/double pthigh = $pt/" applyWeightsFor1GeV_runner.C

	root -l -b -q applyWeightsFor1GeV_runner.C

	# All D0 pT [pt-10]

	cp applyWeightsFor1GeV.C applyWeightsFor1GeV_runner.C
	sed -i'.bak' 's/applyWeightsFor1GeV/applyWeightsFor1GeV_runner/g' applyWeightsFor1GeV_runner.C
	sed -i'.bak' "s/double ptlow = 1/double ptlow = $pt/" applyWeightsFor1GeV_runner.C

	root -l -b -q applyWeightsFor1GeV_runner.C

done