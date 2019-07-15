set term png
set output './dat_images/refs_two_tanks.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [1:3]
plot 'two_tanks_id_fcomp.dat' using 1:2 with lines title 'x1', 'two_tanks_id_fcomp.dat' using 1:3 with lines title 'x2'