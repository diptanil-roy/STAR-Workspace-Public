#!/bin/bash
starver SL16d
make clean
make
./anaEmbeddingEff
make clean
rm AuAu200GeV.KaonMinus.Embedding.Efficiency.pdf
