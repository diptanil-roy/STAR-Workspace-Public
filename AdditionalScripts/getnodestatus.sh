#!/bin/bash

condor_q -global droy1 > /star/u/droy1/condor_output.list
filename="/star/u/droy1/condor_output.list"

printf '%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n' "Node" "Total Jobs" "Active Jobs" "Idle Jobs" "My Total Jobs" "My Active Jobs" "My Idle Jobs" "My Held Jobs"
for node in {6005..6016}
do
	linenum=`grep -n "rcas$node" $filename | awk -F: '{print $1}'`
	totaljobs=`tail -n +$linenum $filename | grep "Total for all users" | head -1 | awk '{print $5}'`
	runningjobs=`tail -n +$linenum $filename | grep "Total for all users" | head -1 | awk '{print $13}'`
	idlejobs=`tail -n +$linenum $filename | grep "Total for all users" | head -1 | awk '{print $11}'`

	mytotaljobs=`tail -n +$linenum $filename | grep "Total for query" | head -1 | awk '{print $4}'`
	myrunningjobs=`tail -n +$linenum $filename | grep "Total for query" | head -1 | awk '{print $12}'`
	myidlejobs=`tail -n +$linenum $filename | grep "Total for query" | head -1 | awk '{print $10}'`
	myheldjobs=`tail -n +$linenum $filename | grep "Total for query" | head -1 | awk '{print $14}'`
	
	printf '%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n' "rcas$node" "$totaljobs" "$runningjobs" "$idlejobs" "$mytotaljobs" "$myrunningjobs" "$myidlejobs" "$myheldjobs" 
done

rm $filename