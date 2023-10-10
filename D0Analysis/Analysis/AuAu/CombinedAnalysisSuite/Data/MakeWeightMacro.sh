#!/bin/bash

[[ -z "$1" ]] && { echo "Enter The Name For the Plots Directory As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Name For the sWeights As Argument #2" ; exit 1; }

mkdir -p $1
mkdir -p $2

for pt in {1..5}
do
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 0, 10)" #Central
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 10, 20)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 20, 30)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 30, 40)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 40, 80)" #Peripheral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", $pt, 10, 0, 80)" #Total
done


lowptbin=(1 2 3 4 5)
highptbin=(2 3 4 5 10)

arraylength=${#lowptbin[@]}

for (( i=0; i<${arraylength}; i++ ));
do
    # echo "${lowptbin[$i]}, ${highptbin[$i]}"
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 0, 10)" #Central
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 10, 20)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 20, 30)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 30, 40)" #MidCentral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 40, 80)" #Peripheral
    root -l -b -q "makeweights.C(\"$1\", \"$2\", ${lowptbin[$i]}, ${highptbin[$i]}, 0, 80)" #Total
done