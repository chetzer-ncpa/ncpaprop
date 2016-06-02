#!/bin/bash

cd ../samples
../bin/pape  --ncpatoy --azimuth 90 --freq 0.1 --write_2D_TLoss

# compare files:
# profile_int.dat
# attn.pe
# tloss_2d.pe
# tloss_1d.pe
mkdir -p ../test/results/pape
mv profile_int.dat attn.pe tloss_2d.pe tloss_1d.pe ../test/results/pape/

