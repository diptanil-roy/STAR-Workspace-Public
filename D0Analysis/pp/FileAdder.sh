#!/bin/bash

#rm *.list

mkdir -p ./AddedFiles
mkdir -p ./log

touch filelist1.list

ls -d ./*.root > filelist1.list

touch filelist2.list

touch difffilelist.list

comm -13 <(sort filelist2.list) <(sort filelist1.list) > difffilelist.list

n=$(wc -l difffilelist.list | awk '{print $1}')

if [[ $n -gt 1000 ]]
then
	timestamp=$(date +%s)
	echo $(date +%H-%M-%S) 
	hadd AddedFiles/$timestamp.root @difffilelist.list > log/$timestamp.log
	cp filelist1.list filelist2.list
else
	timestamp=$(date +%s)
	newestfile=$(ls -t AddedFiles/ | head -n1 | sed 's@AddedFiles/@@' | sed 's@.root@@')
	elapsed=$((timestamp-newestfile))
	if [[ $elapsed -gt 1800 && $n -gt 0 ]]
	then
		echo $(date +%H-%M-%S)
		hadd AddedFiles/$timestamp.root @difffilelist.list > log/$timestamp.log
		cp filelist1.list filelist2.list
	fi
fi

