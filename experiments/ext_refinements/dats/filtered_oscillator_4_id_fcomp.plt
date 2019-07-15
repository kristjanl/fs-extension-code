set term png
set output './dat_images/refs_filtered_oscillator_4.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [7:17]
plot 'filtered_oscillator_4_id_fcomp.dat' using 1:2 with lines title 'x', 'filtered_oscillator_4_id_fcomp.dat' using 1:3 with lines title 'y', 'filtered_oscillator_4_id_fcomp.dat' using 1:4 with lines title 'x1', 'filtered_oscillator_4_id_fcomp.dat' using 1:5 with lines title 'x2', 'filtered_oscillator_4_id_fcomp.dat' using 1:6 with lines title 'x3', 'filtered_oscillator_4_id_fcomp.dat' using 1:7 with lines title 'z'