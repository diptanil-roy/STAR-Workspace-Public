#!/bin/bash

D0pTLow=1

[[ -z "$1" ]] && { echo "Enter The Directory You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Superiterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }
[[ -z "$4" ]] && { echo "Choose what distribution you want to use. Choices now are 0 = FONLL, 1 = FONLL Reweighed With Data" ; exit 1; }

usedatadist=$4

D0pTLow=${5:-1}
echo "D0 pT Low is $D0pTLow"

# rm -rf $1$3

mkdir -p $1$3
mkdir -p $1$3/Central
mkdir -p $1$3/MidCentral
mkdir -p $1$3/Peripheral
mkdir -p $1$3/Plots
mkdir -p $1$3/Plots_RandomizedResponse

## Doing the MC Bit

root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 0, false, $usedatadist)" #Response Matrix Histograms Filled

if [ "$usedatadist" -eq 1 ]
then
	root -l -b -q "CompareDataVsGeant.C(\"$1$3\", $D0pTLow, \"$1$3\", \"After_Reweight\")"
fi

#This step unfolds the measured histograms using the unweighted response matrix
#I had this step for test only, to make sure that my weighing procedure works. The idea is unfolding this vs unfolding after the response matrix goes through the whole process, only with the weights set to 1.
#On April 4, 2023, I found that the unfolding is identical before and after the weighing process for weights = 1. So, I commented this step out.
# root -l -b -q "UnfoldData.C(\"$1$3\", 0, $3, kTRUE)" 

for (( siter = 0; siter <= $2; siter++ ));
do
	root -l -b -q "CreateResponseMatrix.C(\"$1$3\", $siter, $3, 0)"
	if [ "$usedatadist" -eq 6 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, 1)"
	elif [ "$usedatadist" -eq 7 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, -1)"
	elif [ "$usedatadist" -eq 8 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, 2)"
	elif [ "$usedatadist" -eq 9 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, -2)"
	elif [ "$usedatadist" -eq 10 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, 3)"
	elif [ "$usedatadist" -eq 11 ]
	then
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, -3)"
	else
		root -l -b -q "UnfoldData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE)"
	fi
	# root -l -b -q "RandomizeResponse.C(\"$1$3\", $siter, $3, 0)"
	# root -l -b -q "UnfoldUsingRandomData.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE)"
done

/Applications/cpdf `ls -tr $1$3/Central/Step*.pdf` -o $1$3/Central/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/MidCentral/Step*.pdf` -o $1$3/MidCentral/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/Peripheral/Step*.pdf` -o $1$3/Peripheral/SISteps.pdf

ClosureDirectory=${1/Data*/Closure}$3

root -l -b -q "PlotTheUnfoldedDistributionForData.C(\"$1$3\", $D0pTLow, $2, $3, \"$ClosureDirectory/Plots_CS1\", kTRUE)"
root -l -b -q "PlotTheUnfoldedDistributionForAreaData.C(\"$1$3\", $D0pTLow, $2, $3, \"$ClosureDirectory/Plots_AreaBased\", kTRUE)"
# root -l -b -q "PlotTheUnfoldedDistributionForRandomizedResponse.C(\"$1$3\", $2, $3, \"\", kTRUE)"
