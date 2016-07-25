#!/bin/bash

pape --g2senvfile ../../rdmodess_gp_files/code_outputs/g2sgcp2011012606L.jordan.env \
                --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 \
                --freq 0.3 --sourceheight_km 0 --receiverheight_km 0 \
                --maxheight_km 180 --starter_type gaussian --n_pade 6 \
                --maxrange_km 500 --write_2D_TLoss

mv tloss_1d.pe ex1_tloss_1d.pe
mv tloss_2d.pe ex1_tloss_2d.pe 


pape --use_1D_profiles_from_dir ../../rdmodess_gp_files/code_outputs/profiles \
                --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 \
                --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 180 \
                --starter_type gaussian --n_pade 6 --maxrange_km 1000 \
                --use_profile_ranges_km 0_20_60_400 --write_2D_TLoss

mv tloss_1d.pe ex2_tloss_1d.pe
mv tloss_2d.pe ex2_tloss_2d.pe 
