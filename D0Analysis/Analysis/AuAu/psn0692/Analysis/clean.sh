#!/bin/bash
cd vtxCorr  # vtx correction
rm -rf pic data

cd ../default  # raw yield from gaus+pol1 fit method
rm -rf pic data

cd ../count  # raw yield from side-band method
rm -rf pic data

cd ../fitRange # raw yield from gaus+pol1 fit method, fit range was changed
rm -rf pic data

cd ../likeSign # raw yield from sameEvent likeSign background subtraction
rm -rf pic data

cd ../ptCut1  # pT>0.6 ==> pT>0.3
rm -rf pic data

cd ../ptCut2  # pT>0.6 ==> pT>0.5
rm -rf pic data

cd ../topoCut1  # tight topo cuts
rm -rf pic data

cd ../topoCut2  # loose topo cuts
rm -rf pic data

cd ../DoubleCount # double count
rm -rf data

cd ../sys     # sys error summary
rm -rf data pic

cd ../ptShift  # do pT shift
rm -rf pic

cd ../RAA/yifeiFit_USED
rm -rf pic
cd ..

cd ../Rcp
rm -rf pic

cd ../Xsection
rm -rf pic data
