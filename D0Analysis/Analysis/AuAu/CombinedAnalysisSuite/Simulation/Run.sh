#!/bin/bash

D0pTLow=1

[[ -z "$1" ]] && { echo "Enter The Base Directory Name You Want to Put the Files In As Argument #1" ; exit 1; }
[[ -z "$2" ]] && { echo "Enter The Number Of Superiterations as Argument #2" ; exit 1; }
[[ -z "$3" ]] && { echo "Enter The Number Of Iterations as Argument #3" ; exit 1; }
[[ -z "$4" ]] && { echo "Choose what distribution you want to use. Choices now are 0 = FONLL, 1 = FONLL Reweighed With Data" ; exit 1; }
[[ -z "$5" ]] && { echo "No Extra Tag Used" ; }

ExtraTag=$5
# root -l -b -q "MakeFONLLHistogram.C(1, 30)"
# root -l -b -q "MakePYTHIAHistogram.C(1, 30)"

for pt in {1..5}
do
    D0pTLow=$pt
    echo "D0 pT Low is $D0pTLow"

    ./RunFullChain.sh $1 $2 $3 $4 $D0pTLow $ExtraTag &

done

wait

echo "We are done!"