unset key
set style data lines

dir="~/work/atmospheric/infrasound/papers/stratospheric/buncefield_2/figures_workspace/strat_jet_alone/60_ms_jet/"

set parametric

set xlabel "Soundspeed [m/s]
set ylabel "Altitude [km]" offset 1

set arrow from 208,94.5 to 453,94.5 heads lt 3 lw 3
set arrow from 314,10 to 453,10 heads lt 1 lw 3
set title "Relevant Phase Speeds"

set term post enh eps color solid 24
set out "../wvnums_modess.eps"
plot [0:20] [150:500] [0:160] dir."toymodel_60_ms_jet.dat" using (20*sqrt($5)+$2):1 lt 7 lw 5,\
     206.4,t+84.5 lt 3 lw 5,454,t+84.5 lt 3 lw 5,\
     314,t lt 1 lw 5,454,t lt 1 lw 5

unset arrow

! epstopdf ../wvnums_modess.eps;rm ../wvnums_modess.eps; mv wvnums_modess.pdf ../wvnums_modess.pdf