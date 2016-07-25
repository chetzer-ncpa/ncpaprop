unset key 

set pm3d map
set palette defined (0 "white",1 "yellow",2 "red",3 "blue")

unset grid

c(x,y)=x*cos(pi*y/180.0)
s(x,y)=x*sin(pi*y/180.0)

set xlabel "Easterly [km]" 
set ylabel "Northerly [km]"
set xtics rotate
set cblabel "Tloss [dB]" offset 1

set xrange [-1000:1000]
set yrange [-1000:1000]
set cbrange [-130:-90]

set size 0.8,0.8
set term png giant
set size square

set out "../modess_ex3.png"
set title "Transmission Loss Magnitude; 0.1 Hz"

splot "../example_outputs/Nby2D_tloss_1d.nm" \
      using (s($1,$2)):(c($1,$2)):(TL(mag($3,$4)))

set out "../modess_ex3_inco.png"
set title "Incoherent Transmission Loss; 0.1 Hz"

splot "../example_outputs/Nby2D_tloss_1d.nm" \
      using (s($1,$2)):(c($1,$2)):(TL($5))
      
