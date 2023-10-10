#!/bin/bash

#JEWEL to Root
#This script will convert the output of JEWEL to a root file

root -l -b -q 'JEWELToRoot.C("'$1'")'
