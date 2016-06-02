#!/bin/bash

cd ../samples
../bin/raytrace.2d --azimuth 90 --elev 1 --delev 1 --maxelev 45 --skips 1 --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --maxraylength 800 --maxheight 140
cat raypath_az090* > raypaths.2d.dat
rm raypath_az090*

# compare the following files:
# raypaths.2d.dat
mkdir -p ../test/results/raytrace.2d/
mv raypaths.2d.dat ../test/results/raytrace.2d/

