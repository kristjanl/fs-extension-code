set term png
set output './dat_images/refs_filtered_oscillator_8.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [7:17]
plot 'filtered_oscillator_8_id_fcomp.dat' using 1:2 with lines title 'x', 'filtered_oscillator_8_id_fcomp.dat' using 1:3 with lines title 'y', 'filtered_oscillator_8_id_fcomp.dat' using 1:4 with lines title 'f4a_x1', 'filtered_oscillator_8_id_fcomp.dat' using 1:5 with lines title 'f4a_x2', 'filtered_oscillator_8_id_fcomp.dat' using 1:6 with lines title 'f4a_x3', 'filtered_oscillator_8_id_fcomp.dat' using 1:7 with lines title 'f8_x1', 'filtered_oscillator_8_id_fcomp.dat' using 1:8 with lines title 'f4b_x1', 'filtered_oscillator_8_id_fcomp.dat' using 1:9 with lines title 'f4b_x2', 'filtered_oscillator_8_id_fcomp.dat' using 1:10 with lines title 'f4b_x3', 'filtered_oscillator_8_id_fcomp.dat' using 1:11 with lines title 'z'