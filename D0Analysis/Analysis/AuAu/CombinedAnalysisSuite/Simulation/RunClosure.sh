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
mkdir -p $1$3/Plots_AreaBased_DataUnc
mkdir -p $1$3/Plots_AreaBased
mkdir -p $1$3/Plots_CS1_DataUnc
mkdir -p $1$3/Plots_CS1

# Generating the Data Dependent Weights

# # Doing the MC Bit

root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 1, false, $usedatadist)" #Split 1

root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 2, false, $usedatadist)" #Split 2

if [ "$usedatadist" -eq 1 ]
then
	root -l -b -q "CompareDataVsGeant.C(\"$1$3\", $D0pTLow, \"$1$3\", \"AfterReweight_Split1\")"

	root -l -b -q "CompareDataVsGeant.C(\"$1$3\", $D0pTLow, \"$1$3\", \"AfterReweight_Split2\")"
fi

root -l -b -q "CreateResponseMatrix.C(\"$1$3\", 0, $3, 2)" #Split 2
root -l -b -q "CreateResponseMatrix.C(\"$1$3\", 1, $3, 2)" #Split 2

for (( siter = 0; siter <= $2; siter++ ));
do
	# root -l -b -q "UnfoldMC.C(\"$1$3\", $D0pTLow, $siter, $3, kFALSE, kFALSE)"
	root -l -b -q "UnfoldMC.C(\"$1$3\", $D0pTLow, $siter, $3, kFALSE, kTRUE)"

	# Mimicking Data Uncertainty Here
	#root -l -b -q "UnfoldMC.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, kFALSE)"
	#root -l -b -q "UnfoldMC.C(\"$1$3\", $D0pTLow, $siter, $3, kTRUE, kTRUE)"
done

/Applications/cpdf `ls -tr $1$3/Central/Step*.pdf` -o $1$3/Central/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/MidCentral/Step*.pdf` -o $1$3/MidCentral/SISteps.pdf
/Applications/cpdf `ls -tr $1$3/Peripheral/Step*.pdf` -o $1$3/Peripheral/SISteps.pdf

# root -l -b -q "PlotTheUnfoldedDistributionForMC.C(\"$1$3\", $D0pTLow, $2, $3, kFALSE, kFALSE)"
root -l -b -q "PlotTheUnfoldedDistributionForMC.C(\"$1$3\", $D0pTLow, $2, $3, kFALSE, kTRUE)"
# root -l -b -q "PlotTheUnfoldedDistributionForMC.C(\"$1$3\", $D0pTLow, $2, $3, kTRUE, kFALSE)"
# root -l -b -q "PlotTheUnfoldedDistributionForMC.C(\"$1$3\", $D0pTLow, $2, $3, kTRUE, kTRUE)"