#!/bin/bash

#root -l -b -q 'UnfoldPtvZ.C(1)'
# root -l -b    'UnfoldPtvZ.C(2)'

# root -l -q 'ResponseMatrix_Weighted.C(200.0)' #200 is weighing with the peripheral case only
# root -l -q 'ResponseMatrix_Weighted.C(300.0)' #300 is weighing each centrality case with its own
# # root -l -q 'ResponseMatrix_Weighted.C(2.0)'
# # root -l -q 'ResponseMatrix_Weighted.C(3.0)'
# # root -l -q 'ResponseMatrix_Weighted.C(4.0)'
# # root -l -q 'ResponseMatrix_Weighted.C(5.0)'

root -l -b -q 'ClosureMacro.C(1)'
root -l -b    'ClosureMacro.C(2)'