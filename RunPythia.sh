#!/bin/bash

#WriteXSecWeightsFile = $1

export PtHatBinLows=(15 20 25 35 45 55 -1)

export doWeakDecays="off"

export WriteXSecWeightsFile="XSecGenWeightspSet17Gluons.dat"

rm "$WriteXSecWeightsFile"

touch "$WriteXSecWeightsFile"

for pthbin in 0 1 2 3 4 5
do
	echo "In PtHat bin ${PtHatBinLows[$pthbin]} to ${PtHatBinLows[$((pthbin+1))]}"
	root -b -q -l "PythiaTest.C(${PtHatBinLows[$pthbin]}, ${PtHatBinLows[$((pthbin+1))]}, \"$doWeakDecays\", 0, \"$WriteXSecWeightsFile\")"
done
