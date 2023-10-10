#!/bin/bash

D0pTLow=1

[[ -z "$1" ]] && { echo "Enter The Directory You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Superiterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }

D0pTLow=${4:-1}
echo "D0 pT Low is $D0pTLow"

rm -rf $1$3
# 
mkdir -p $1$3
mkdir -p $1$3/Central
mkdir -p $1$3/MidCentral
mkdir -p $1$3/Peripheral
mkdir -p $1$3/Plots
# 
## Doing the MC Bit
# 
root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 1)" #Split 1
root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 2)" #Split 2
# 
root -l -b -q "CreateResponseMatrix.C(\"$1$3\", 0, $3, 1)" #Split 1
root -l -b -q "CreateResponseMatrixWithDataZ.C(\"$1$3\", $D0pTLow, 0, $3, 2)" #Split 2
root -l -b -q "CreateResponseMatrixWithDataZ.C(\"$1$3\", $D0pTLow, 1, $3, 2)" #Split 2
# 
for (( siter = 0; siter <= $2; siter++ ));
do
	root -l -b -q "CreateResponseMatrix.C(\"$1$3\", $siter, $3, 1)" #Split 1 - Flattened Z
	root -l -b -q "UnfoldMC.C(\"$1$3\", $siter, $3, kFALSE)"
done
# 
/Applications/cpdf `ls -tr $1$3/Central/Step*.pdf` -o $1$3/Central/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/MidCentral/Step*.pdf` -o $1$3/MidCentral/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/Peripheral/Step*.pdf` -o $1$3/Peripheral/SISteps.pdf

root -l "PlotTheUnfoldedDistributionForMC.C(\"$1$3\", $2, $3)"