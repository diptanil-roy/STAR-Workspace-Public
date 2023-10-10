#!/bin/bash

for i in {3..7}
do
   root -l -b -q "HIOverlayMacro.C($i)" &
   sleep 60
done

# wait

# for i in {12..20}
# do
#    root -l -b -q "HIOverlayMacro.C($i)" &
#    sleep 60
# done

