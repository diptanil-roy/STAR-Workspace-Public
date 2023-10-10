#!/bin/bash

[[ -z "$1" ]] && { echo "Enter The Base Directory Name You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Superiterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }
[[ -z "$4" ]] && { echo "Choose what distribution you want to use. Choices now are 0 = FONLL, 1 = FONLL Reweighed With Data" ; exit 1; }

usedatadist=$4

D0pTLow=$5
echo "D0 pT Low is $D0pTLow"

if [[ -z "$6" ]]
then
   ExtraTag=""
else
   ExtraTag="_"$6
fi

echo "Extra Tag is $ExtraTag"

# root -l -b -q "MakeFONLLHistogram.C($D0pTLow, 30)"
# root -l -b -q "MakePYTHIAHistogram.C($D0pTLow, 30)"

CentDir="CentWeight_D0pT_$D0pTLow"GeV
rm -rf $CentDir
mkdir -p $CentDir

ClosureDirectory=$1"_"$D0pTLow"GeV_Closure"$ExtraTag
DataDirectory=$1"_"$D0pTLow"GeV_Data"$ExtraTag

BaseDirectory=$1"_"$D0pTLow"GeV_Data"

rm -rf $ClosureDirectory$3
rm -rf $DataDirectory$3

mkdir -p $ClosureDirectory$3
mkdir -p $DataDirectory$3

if [ "$usedatadist" -eq 1 ]
then
	root -l -b -q "CompareDataVsGeant.C(\"$BaseDirectory$3\", $D0pTLow, \"$CentDir\")"
   root -l -b -q "CompareDataVsGeant.C(\"$BaseDirectory$3\", $D0pTLow, \"$ClosureDirectory$3\")"
   root -l -b -q "CompareDataVsGeant.C(\"$BaseDirectory$3\", $D0pTLow, \"$DataDirectory$3\")"
fi

root -l -b -q "PrepareResponseHistograms.C(\"$CentDir\", $D0pTLow, 0, true, $usedatadist)"
root -l -b -q "GetCentralityWeights.C($D0pTLow)"

./RunClosure.sh $ClosureDirectory $2 $3 $usedatadist $D0pTLow

./RunDataUnfolding.sh $DataDirectory $2 $3 $usedatadist $D0pTLow

echo "Done!"
echo "Closure Directory is $ClosureDirectory$3"
echo "Data Directory is $DataDirectory$3"
