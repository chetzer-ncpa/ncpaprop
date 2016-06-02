#!/bin/bash

cd ../samples
../bin/CModess --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1

# compare the following files:
# tloss_1d.cnm
# tloss_1d.lossless.cnm
mkdir -p ../test/results/CModess
mv tloss_1d.cnm tloss_1d.lossless.cnm ../test/results/CModess
