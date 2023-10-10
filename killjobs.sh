#!/bin/bash

while read -r line; do condor_rm $line; done < jobstokill
