#!/bin/bash

[[ -z "$1" ]] && { echo "Enter The Directory You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Iterations as Argument #2" ; exit 1; }

mkdir -p $1$2
mkdir -p $1$2/Central
mkdir -p $1$2/MidCentral
mkdir -p $1$2/Peripheral

root -l -b -q "CreateOldRespMacro.C(1, $2, \"$1$2\", 1)"
root -l -b -q "CreateOldRespMacro.C(1, $2, \"$1$2\", 2)"
root -l -b -q "UnfoldSimOldWayMacro.C(1, $2, \"$1$2\", \"true\")"

root -l -b -q "CreateOldRespMacro.C(1, $2, \"$1$2\")"
root -l -b -q "UnfoldDataOldWayMacro.C(1, $2, \"$1$2\", \"true\", \"Resp1D\")"
root -l -b -q "MakeOldDataComparisonPlots.C($2, \"$1$2\", \"true\", 1)"