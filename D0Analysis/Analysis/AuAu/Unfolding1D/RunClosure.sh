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

# Doing the MC Bit

root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 1)" #Split 1
root -l -b -q "PrepareResponseHistograms.C(\"$1$3\", $D0pTLow, 2)" #Split 2

root -l -b -q "UnfoldMC.C(\"$1$3\", 0, $3, kFALSE)"
