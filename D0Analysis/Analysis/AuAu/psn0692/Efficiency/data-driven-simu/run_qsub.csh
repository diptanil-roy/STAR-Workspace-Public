#!/bin/tcsh
# set nevents=${nevents1}
# set particle=${particle1}
# set runIndex=${runIndex1}
set nevents=$1
set particle=$2
set runIndex=$3

cd /global/homes/x/xgn1992/rnc_Global/SL16dRun14/myAnalysis_QM17_Jan17/Xiaolong/Redo_0721/efficiency/simu_combine_extendpT
starver SL16d  # this star version should keep same with you complie by ".L toyMcBtoD.C++
root4star -b -q runBtoD.C\(${nevents},\"\/global\/homes\/x\/xgn1992\/rnc_Global\/SL16dRun14\/myAnalysis_QM17_Jan17\/Xiaolong\/Redo_0721\/efficiency\/simu_combine_extendpT\/out\/${particle}_${runIndex}\",\"${particle}\"\)
