#!/bin/bash

[[ -z "$1" ]] && { echo "Enter The Directory You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of SuperIterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }

mkdir -p $1$3
mkdir -p $1$3/Central
mkdir -p $1$3/MidCentral
mkdir -p $1$3/Peripheral

for (( siter=0; siter<=$2; siter++ ))
do
	echo "SUPERITER $siter"

	root -l -b -q "PrepareResponseMatrix.C($siter, $3, \"$1$3\", false, true)"

	root -l -b -q "CreateRespMacro.C($siter, $3, \"$1$3\")"

	root -l -b -q "UnfoldDataMacro.C($siter, $3, \"$1$3\")"
done

root -l -b -q "MakeDataComparisonPlots.C($3, \"$1$3\", $2)"

convert -layers OptimizePlus -delay 50 $(ls $1$3/Central/*.png | sort -V) $1$3/Central/Central.gif
convert -layers OptimizePlus -delay 50 $(ls $1$3/MidCentral/*.png | sort -V) $1$3/MidCentral/MidCentral.gif
convert -layers OptimizePlus -delay 50 $(ls $1$3/Peripheral/*.png | sort -V) $1$3/Peripheral/Peripheral.gif