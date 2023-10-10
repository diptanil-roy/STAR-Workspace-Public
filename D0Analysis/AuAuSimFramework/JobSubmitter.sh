#!/bin/bash

mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/HIOverlay_HFJets_Apr12_2023_NoJetpTCutOff_pthat_3_inf
mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/HIOverlay_HFJets_Apr12_2023_NoJetpTCutOff_pthat_3_inf/gen
mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/HIOverlay_HFJets_Apr12_2023_NoJetpTCutOff_pthat_3_inf/log
mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/HIOverlay_HFJets_Apr12_2023_NoJetpTCutOff_pthat_3_inf/err
mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/HIOverlay_HFJets_Apr12_2023_NoJetpTCutOff_pthat_3_inf/out

sort -R Run14MidLumi_Randomised_2.list > RandomizedRun14List.list

for i in {3..200}
do
	echo "Round $i"

	a=$(( 100*($i-1) + 1 ))

	tail -n +$a RandomizedRun14List.list | head -n +100 > $i.list

	cp HIResponse_3_inf_Full.xml $i.xml

	sed -i "s/RandomizedRun14List.list/$i.list/" $i.xml

	star-submit $i.xml

	rm $i.list
	rm $i.xml
done

