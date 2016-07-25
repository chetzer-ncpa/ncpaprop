#!/bin/bash

ModessRD1WCM --use_1D_profiles_from_dir profiles \
             --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 \
             --use_profile_ranges_km 50_100_150_200_250 --write_2D_TLoss

mv tloss_rd_1d.lossless.nm ex1_tloss_rd_1d.lossless.nm
mv tloss_rd_1d.nm ex1_tloss_rd_1d.nm 
mv tloss_rd_2d.nm ex1_tloss_rd_2d.nm 

ModessRD1WCM --g2senvfile g2sgcp2011012606L.jordan.env \
             --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 \
             --use_profile_ranges_km 50_100_150_200_250 --write_2D_TLoss

mv tloss_rd_1d.lossless.nm ex2_tloss_rd_1d.lossless.nm
mv tloss_rd_1d.nm ex2_tloss_rd_1d.nm 
mv tloss_rd_2d.nm ex2_tloss_rd_2d.nm 

