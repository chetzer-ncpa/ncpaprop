set key samplen 2

dir="../example_outputs/"

set xlabel "Range [km]"
set ylabel "TL [dB]"

set title "1D Transmission Loss Magnitude; 0.1 Hz, East"
set term post enh eps color solid 22
set out "cmodess_ex1.eps"
plot dir."tloss_1d.cnm" using 1:(TL(mag($2,$3))) lt 3 lw 3 title "CModess",\
     dir."tloss_1d.cnm" using 1:(TL($4)) lt 2 lw 4 title "CModess inco"

! epstopdf cmodess_ex1.eps;rm cmodess_ex1.eps;mv cmodess_ex1.pdf ..


set title "1D Transmission Loss Magnitude; 0.1 Hz, North"
set term post enh eps color solid 22
set out "1D_cmodess_ex_N.eps"
plot dir."tloss_1d_N.cnm" using 1:(TL(mag($2,$3))) lt 3 lw 3 title "CModess",\
     dir."tloss_1d_N.cnm" using 1:(TL($4)) lt 2 lw 4 title "CModess inco"

! epstopdf 1D_cmodess_ex_N.eps;rm 1D_cmodess_ex_N.eps;mv 1D_cmodess_ex_N.pdf ..
