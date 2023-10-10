#!/bin/bash

ln -s ../DataAna/D0_Data_Mix_MergeUsed.root .

cd vtxCorr  # vtx correction
rm -rf pic data
bash run.sh

cd ../default  # raw yield from gaus+pol1 fit method
rm -rf pic data
bash run.sh

cd ../count  # raw yield from side-band method
rm -rf pic data
bash run.sh

cd ../fitRange # raw yield from gaus+pol1 fit method, fit range was changed
rm -rf pic data
bash run.sh

cd ../likeSign # raw yield from sameEvent likeSign background subtraction
rm -rf pic data
bash run.sh

cd ../default
root -b -q re_get_yield.C # re calculate the yield for some pt/centrality bin
root -b -q get_Rcp.C

cd ../count
root -b -q get_sys.C

cd ../fitRange
root -b -q get_sys.C

cd ../likeSign
root -b -q get_sys.C

cd ../ptCut1  # pT>0.6 ==> pT>0.3
rm -rf pic data
bash run.sh
root -b -q get_Rcp.C

cd ../ptCut2  # pT>0.6 ==> pT>0.5
rm -rf pic data
bash run.sh
root -b -q get_Rcp.C

cd ../topoCut1  # tight topo cuts
rm -rf pic data
bash run.sh
root -b -q get_Rcp.C

cd ../topoCut2  # loose topo cuts
rm -rf pic data
bash run.sh
root -b -q get_Rcp.C

cd ../DoubleCount # double count
rm -rf data
root -b -q rm_doubleCount.C

cd ../sys     # sys error summary
rm -rf data pic
root -b -q write_err.C
root -b -q plot_err.C
root -b -q write_Rcp_err.C
root -b -q plot_Rcp_err.C
root -b -q write_CrossS_err.C

cd ../ptShift  # do pT shift
rm -rf pic
root -b -q doPtShift.C
root -b -q plot_spectra.C

# cd ../RAA
cd ../RAA/yifeiFit_USED
rm -rf pic
# root -b -q write_RAA_Integral.C # use RAA in integral pT bin,  write_RAA.C ==> use RAA at pT value after pT shift
root -b -q write_RAA.C # use RAA in integral pT bin,  write_RAA.C ==> use RAA at pT value after pT shift
root -b -q plot_RAA_1.C  # draw each cent
root -b -q plot_RAA.C    # draw 0-10%, 10-40%, 40-80% together
root -b -q getfitsys_xx.C
root -b -q write_RAA.C # use RAA in integral pT bin,  write_RAA.C ==> use RAA at pT value after pT shift
cd ..

cd ../Rcp
rm -rf pic
root -b -q plot_Rcp.C
root -b -q plot_Rcp_fit.C
root -b -q plot_Rcp1_pTshift.C
root -b -q plot_Rcp2_pTshift.C

cd ../Xsection
root -b -q plot_D0CrossS.C
root -b -q plot_D0CrossS2.C
root -b -q plot_D0CrossS_forRcp.C
root -b -q plot_D0CrossS2_forRcp.C
