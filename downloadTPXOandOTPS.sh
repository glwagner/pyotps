#!/bin/bash

# This script downloads and compiles the Oregon Tidal Prediction
# Software (OTPS) and then downloads the TPXO8 data the drives it.

ftpaddr="ftp://ftp.oce.orst.edu"

otpsdir="/dist/tides"
tpxodir="/dist/tides/TPXO8_compact"

otpsname="OTPS2"
tpxoname="tpxo8_atlas_compact_v1"

suffix=".tar.Z"

installdir=$PWD

# Download and untar OTPS 
cd $installdir

wget "$ftpaddr$otpsdir/$otpsname$suffix" ./
tar xvf "$otpsname$suffix"

cd $otpsname

read -p "You are about to download the TPXO global solution, 
which packed in a 1.24 Gigabyte tar file! Press enter to continue."

wget "$ftpaddr$tpxodir/$tpxoname$suffix"
tar xvf "$tpxoname$suffix"

echo "Compiling the extract_HC and predict_tide routines
provided with OTPS. The compilation will fail if there is no 
working version of gfortran in the path."

make extract_HC
make predict_tide
