#myTDpape_dir="/home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic"
myTDpape_dir="testrun"

#./tdpape1 --pulse_prop_src2rcv $myTDpape_dir --range_R_km 240 --waveform_out_file mywavef1.dat --max_celerity 320 --use_builtin_pulse

./tdpape --pulse_prop_src2rcv_grid  $myTDpape_dir  --R_start_km 200 --R_end_km 240 --DR_km 20 --waveform_out_file wfgrid4.dat --max_celerity 320 --src_waveform_file source_waveform.dat
