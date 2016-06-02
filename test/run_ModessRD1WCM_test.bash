#!/bin/bash

cd ../samples
../bin/ModessRD1WCM --use_1D_profiles_from_dir profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profile_ranges_km 100_200

# Compare the following files:
# tloss_rd_1d.lossless.nm
# tloss_rd_1d.nm
mkdir -p ../test/results/ModessRD1WCM
mv tloss_rd_1d.lossless.nm tloss_rd_1d.nm ../test/results/ModessRD1WCM/
