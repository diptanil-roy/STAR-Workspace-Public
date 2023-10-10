#!/bin/bash

touch MissedFilesList_Mar24.list

while read p; do
	less /gpfs01/star/pwg_tasks/jetcorr03/HFJets_Mar24_2022_SimulatedJets/gen/$p >> MissedFilesList_Mar24.list 
done < UnfinishedJobs_NewSample_Mar24.log

