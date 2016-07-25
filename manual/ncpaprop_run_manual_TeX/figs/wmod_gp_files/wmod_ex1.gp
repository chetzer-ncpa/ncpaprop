set key samplen 2

dir="../example_outputs/"

set xlabel "Range [km]"
set ylabel "TL [dB]"

set title "1D Transmission Loss Magnitude; 0.1 Hz"
set term post enh eps color solid 22
set out "wmod_ex1.eps"
plot dir."wtloss_1d.lossless.nm" using 1:(TL(mag($2,$3))) lt 3 lw 3 title "lossless",\
     dir."wtloss_1d.nm" using 1:(TL(mag($2,$3))) lt 1 lw 3 title "lossy",\
     dir."wtloss_1d.nm" using 1:(TL($4)) lt 2 lw 4 title "lossy inco"

! epstopdf wmod_ex1.eps;rm wmod_ex1.eps;mv wmod_ex1.pdf ..
