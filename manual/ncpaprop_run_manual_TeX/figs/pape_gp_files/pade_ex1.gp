set key samplen 2

dir="code_outputs/"

set xlabel "Range [km]"
set ylabel "TL [dB]"

set title "1D Transmission Loss Magnitude; 0.1 Hz"
set term post enh eps color solid 22
set out "pade_ex1_1d.eps"
plot dir."ex1_tloss_1d.pe" using 1:(TL(mag($2,$3))) lt 3 lw 3 title ""

! epstopdf pade_ex1_1d.eps;rm pade_ex1_1d.eps
! mv pade_ex1_1d.pdf ..
! cp ../pade_ex1_1d.pdf ../new_figs/pade_ex1_1d.pdf


unset key 

set pm3d map
set palette defined (0 "white",1 "yellow",2 "red",3 "blue")

set xtics 100
set ytics 50
set cbtics 10
set xrange [0:1000]
set yrange [0:150]
set cbrange [-130:-90]
set xlabel "Range [km]"
set ylabel "Altitude [km]"
set cblabel "Tloss [dB]" offset 1

set title "2D Transmission Loss Magnitude; 0.1 Hz"
set term png giant
set size 0.8,0.8
set out "pade_ex1_2d.png"
splot [0:500] "code_outputs/ex1_tloss_2d.pe" using 1:2:(TL(mag($3,$4)))

!mv pade_ex1_2d.png ../pade_ex1_2d.png
!cp ../pade_ex1_2d.png ../new_figs/pade_ex1_2d.png

