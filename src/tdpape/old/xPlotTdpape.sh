
# tdpape output: plot propagated pulse to the screen; 
# the second argument is the waveform filename, the third is the range in km
# e.g. 
# ./xPlotTdpape.sh plotpulse mywavf.dat 240
if 
    test $1 = plotpulse
then
  gnuplot -persist \
    -e "set grid; show grid; set xlabel 'Time [s]'; set ylabel 'Amplitude';" \
    -e "set title 'Propagated pulse to $3 km';" \
    -e "plot './$2' using 1:2 with lines lt 1 title '';"
   
  echo $1 $2 $3
fi

# tdpape output: plot propagated pulse to .png
# the second argument is the waveform filename, 
# the third argument is the range in km
# 4th argument is the output png file
# e.g. 
# ./xPlotTdpape.sh plotpulsepng mywavf.dat 240 a.png 
if 
    test $1 = plotpulsepng
then
  gnuplot -persist \
    -e "set term png size 900,500;" \
    -e "set output '$4';" \
    -e "set grid; show grid; set xlabel 'Time [s]'; set ylabel 'Amplitude';" \
    -e "set title 'Propagated pulse to $3 km';" \
    -e "plot './$2' using 1:2 with lines lt 1 title '';"
   
  echo $1 $2
fi
