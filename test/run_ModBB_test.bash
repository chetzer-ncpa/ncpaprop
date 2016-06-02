#!/bin/bash

cd ../samples
../bin/ModBB --out_disp_src2rcv_file myDispersionFile.dat --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --f_step 0.001953125 --f_max 0.5 --use_modess
../bin/ModBB --out_dispersion_files disprs --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --f_step 0.001953125 --f_max 0.5 --use_modess

# compare the following files:
# myDispersionFile.dat
# various disprs_#####_nm.bin
mkdir -p ../test/results/ModBB
mv myDispersionFile.dat disprs_*_nm.bin ../test/results/ModBB


