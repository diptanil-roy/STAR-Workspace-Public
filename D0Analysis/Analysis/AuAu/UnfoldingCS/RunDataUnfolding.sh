#!/bin/bash

D0pTLow=1

[[ -z "$1" ]] && { echo "Enter The Directory You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Superiterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }

D0pTLow=${4:-1}
echo "D0 pT Low is $D0pTLow"

rm -rf $1$3

mkdir -p $1$3
mkdir -p $1$3/Central
mkdir -p $1$3/MidCentral
mkdir -p $1$3/Peripheral
mkdir -p $1$3/Plots

### Doing the MC Bit

root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 0)" #Response Matrix Histograms Filled

#This step unfolds the measured histograms using the unweighted response matrix
#I had this step for test only, to make sure that my weighing procedure works. The idea is unfolding this vs unfolding after the response matrix goes through the whole process, only with the weights set to 1.
#On April 4, 2023, I found that the unfolding is identical before and after the weighing process for weights = 1. So, I commented this step out.
# root -l -b -q "UnfoldData.C(\"$1$3\", 0, $3, kTRUE)" 

for (( siter = 0; siter <= $2; siter++ ));
do
	# root -l -b -q "CreateResponseMatrixExternalWeights.C(\"$1$3\", $siter, $3, 0, \"/Volumes/WorkDrive/STAR-Workspace/D0Analysis/LocalPythiaSim/MacroForSTAR/BBBar.root\")" #This step unfolds the measured histograms after the weighing process.
	root -l -b -q "CreateResponseMatrix.C(\"$1$3\", $siter, $3, 0)"
	root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE)"
done

/Applications/cpdf `ls -tr $1$3/Central/Step*.pdf` -o $1$3/Central/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/MidCentral/Step*.pdf` -o $1$3/MidCentral/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/Peripheral/Step*.pdf` -o $1$3/Peripheral/SISteps.pdf

root -l "PlotTheUnfoldedDistributionForData.C(\"$1$3\", $2, $3, \"\", kTRUE)"