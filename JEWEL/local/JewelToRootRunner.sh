#!/bin/bash

ls ./JEWEL_Files/*.hepmc > ./Filelist.list

# Run the JEWEL to Root converter

while read -r line
do
    root -l -b -q 'JEWELToRoot.C("'$line'")'
done < ./Filelist.list
