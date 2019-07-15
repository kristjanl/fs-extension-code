set term png
set output './dat_images/refs_rod_reactor.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:5]
plot 'rod_reactor_id_fcomp.dat' using 1:2 with lines title 'x', 'rod_reactor_id_fcomp.dat' using 1:3 with lines title 'c1', 'rod_reactor_id_fcomp.dat' using 1:4 with lines title 'c2'