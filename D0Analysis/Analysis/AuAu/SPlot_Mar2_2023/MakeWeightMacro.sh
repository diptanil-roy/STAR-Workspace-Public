#!/bin/bash

mkdir -p plots_D0JetpT_0_100GeV

### Central

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### MidCentral (10 - 20)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=10;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=20;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### MidCentral (20 - 30)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=20;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=30;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### MidCentral (30 - 40)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=30;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=40;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### Peripheral (40 - 60)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=40;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=60;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### Peripheral (60 - 80)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=60;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=80;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

## Peripheral (40 - 80)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=40;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=80;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### All Centralities (0 - 80)

cp makeWeightsFor1GeV.C makeWeights_runner.C
sed -i'.bak' 's/makeWeightsFor1GeV/makeWeights_runner/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_low=0/int _cen_low=0;/g' makeWeights_runner.C
sed -i'.bak' 's/int _cen_high=10/int _cen_high=80;/g' makeWeights_runner.C

### 1 < D0 pT < 10
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 2
sed -i'.bak' 's/double _pt_low=1./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=10./double _pt_high=2./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 2 < D0 pT < 3
sed -i'.bak' 's/double _pt_low=1./double _pt_low=2./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=2./double _pt_high=3./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 3 < D0 pT < 4
sed -i'.bak' 's/double _pt_low=2./double _pt_low=3./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=3./double _pt_high=4./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 4 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=3./double _pt_low=4./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=4./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 1 < D0 pT < 5
sed -i'.bak' 's/double _pt_low=4./double _pt_low=1./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=5./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

### 5 < D0 pT < 10
sed -i'.bak' 's/double _pt_low=1./double _pt_low=5./g' makeWeights_runner.C
sed -i'.bak' 's/double _pt_high=5./double _pt_high=10./g' makeWeights_runner.C
root -l -b -q makeWeights_runner.C

mkdir -p sWeights_D0JetpT_0_100GeV
mv News*.root sWeights_D0JetpT_0_100GeV/.