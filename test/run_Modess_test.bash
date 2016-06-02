#!/bin/bash

cd ../samples
../bin/Modess --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --azimuth 90 --freq 0.1

# compare the following files:
# att_coeff.nm
# ceff.nm
# tloss_1d.lossless.nm
# tloss_1d.nm
mkdir -p ../test/results/Modess/
mv tloss_1d.nm ceff.nm tloss_1d.lossless.nm ../test/results/Modess/
