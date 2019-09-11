set term png
set output './dat_images/refs_cruise_control.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:2]
plot 'cruise_control_id_fcomp.dat' using 1:2 with lines title 'v', 'cruise_control_id_fcomp.dat' using 1:3 with lines title 't'