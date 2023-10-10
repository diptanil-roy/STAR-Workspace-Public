#!/bin/bash

# ./Run.sh Aug14_FONLL 0 3 0
# wait
# # ./Run.sh Aug14_FONLL 0 4 0
# # wait
# ./Run.sh Aug14_FONLL 0 5 0
# wait
# ./Run.sh Aug14_FONLL 0 3 1 WeighedByData
# wait
# # ./Run.sh Aug14_FONLL 0 4 1 WeighedByData
# # wait
# ./Run.sh Aug14_FONLL 0 5 1 WeighedByData
# wait
# ./Run.sh Aug14_PYTHIA 0 3 3
# wait
# ./Run.sh Aug14_PYTHIA 0 4 3
# wait
# ./Run.sh Aug14_PYTHIA 0 5 3
# wait

# ./Run.sh Aug14_FONLL_D0YieldUpper 0 4 6
# wait
# ./Run.sh Aug14_FONLL_D0YieldLower 0 4 7
# wait

./Run.sh Aug14_FONLL_D0YieldNoVCEffUpper 0 4 8
wait
./Run.sh Aug14_FONLL_D0YieldNoVCEffLower 0 4 9
wait
./Run.sh Aug14_FONLL_D0YieldVCEffUpper 0 4 10
wait
./Run.sh Aug14_FONLL_D0YieldVCEffLower 0 4 11
wait

echo "We are done!"