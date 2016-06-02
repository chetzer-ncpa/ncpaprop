# This bash script is an example on how to compute in batch mode 
# single frequency output PE (pape) files 
# which can be subsequently used with the time-domain PE (tdpape) to obtain
# waveforms propagated to various distances. "tdpape" is thus a Fourier
# synthesizer of single-frequency results. See tdpape --help for details
# on obtaining the time-domain propagated waveforms.
#
# The user should change the "pape" command in the main loop below 
# to fit his/her needs. The example here contains this command:
#
##   ../../bin/pape --atmosfile1d ../NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq $freq --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 150 --starter_type gaussian --n_pade 6 --maxrange_km $maxrange_km &
#
#
# As input the user must provide the anticipated maximum range of propagation,
# the maximum frequency of the source waveform, fmax, and the time duration 
# (i.e. time record length T on which the waveform is to be represented).
# The number of necessary frequencies Nfreq is determined by fmax*T.
# The output single-frequency files containing the complex pressure have a 
# naming convention "file_counter_papeTL_frequency". For example:
# "084_papeTL_0.164062" is the file number 84 obtained for a frequency offset
# 0.164062 Hz. The number of files should obviously be Nfreq as expected by
# tdpape for pulse propagation. 

# ---------------------------------------------------------------------------
# User input (probably a good idea to avoid blanks around the equal sign):
# ---------------------------------------------------------------------------
maxrange_km=500;         # maximum propagation range in km
fmax=0.5;                # maximum frequency in the source spectrum
T=512;                   # time record length in seconds
directory="currPapeBB4"; # directory where the pape output is to be placed

#-----------------------------------------------------------------------------

if [ -d $directory ]; then
  echo "Directory $directory exists"
else
  echo "$directory does not exist; making it ..."
  mkdir $directory
fi


j=$(ls $directory | wc -l)  # check if dir is empty
echo $j

if [ $j -ne 0 ]; then
echo "$j  $directory is not empty; empty it before running this script"
exit 0
fi

# obtain Nfreq
Nfreq=`echo "scale=10; $fmax*$T" | bc`
Nfreq=$(printf %.0f $Nfreq) # round Nfreq
echo Nfreq = $Nfreq
fstep=`echo "scale=10; $fmax / $Nfreq" | bc`
echo fstep = $fstep


start_time=`date`

for i in $(seq 1 $Nfreq); do
#for i in $(seq 1 ); do
  dd=$(printf tempp$i) # temporary directories
  echo $dd
    # echo $i $(printf temp$i)
  if [ -d $dd ]; then
    echo dir $dd exists
  else
    mkdir $dd
  fi
  
  freq=`echo "scale=10; $i*$fstep" | bc`
  cd $dd
  echo freq = $freq
  
  ## note: new line needed in next command for bash to return prompt after launching background processes
  #cmnd=$( printf "%s\n" "../../bin/pape  --ncpatoy --azimuth 90 --freq $freq --maxrange_km $maxrange_km" )
  #$cmnd
  
  ## alternate command
  #../../bin/pape  --ncpatoy --azimuth 90 --freq $freq --maxrange_km $maxrange_km &
  
  # the user should change this command
  ../../../bin/pape --atmosfile1d ../NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq $freq --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 150 --starter_type gaussian --n_pade 6 --maxrange_km $maxrange_km &

  #`echo "$cmnd &"`

  echo iteration $i
   
  cd ../
 
done 

wait


for i in $(seq 1 $Nfreq); do
  dd=$(printf tempp$i) # temporary directories
  echo $dd
    # echo $i $(printf temp$i)
  
  freq=`echo "scale=10; $i*$fstep" | bc`
  cd $dd
  echo freq = $freq
  
  echo got here $i
  fn=$(printf %03d_papeTL_%g $i $freq)
  echo $fn
  mv tloss_1d.pe $fn
  mv $fn ../$directory  # copy file to directory defined at the beginning
   
  cd ../
  rm -r "$dd" # remove temporary directories
  
done 

end_time=`date`
echo start_time: $start_time
echo end_time: $end_time


