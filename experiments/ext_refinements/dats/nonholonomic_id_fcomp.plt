set term png
set output './dat_images/refs_nonholonomic.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:2]
plot 'nonholonomic_id_fcomp.dat' using 1:2 with lines title 'x', 'nonholonomic_id_fcomp.dat' using 1:3 with lines title 'y', 'nonholonomic_id_fcomp.dat' using 1:4 with lines title 'z'