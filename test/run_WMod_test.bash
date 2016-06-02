#!/bin/bash

cd ../samples
../bin/WMod --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1


# compare the following files:
# wtloss_1d.lossless.nm
# wtloss_1d.nm
mkdir -p ../test/results/WMod/
mv wtloss_1d.lossless.nm wtloss_1d.nm ../test/results/WMod/
