#!/bin/bash

touch tmp.list
grep -l "done with output/error channels" log/*.log >> tmp.list

touch missedfileslist.list

while read p; do
	filename=$(echo $p | sed 's/.log/.list/' | sed 's@log/@gen\/ppJetTreeJob@')
	echo $filename
	cat $filename >> missedfileslist.list
done < tmp.list

rm tmp.list