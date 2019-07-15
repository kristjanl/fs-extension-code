set term png
set output './dat_images/refs_and_v3.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:11]
plot 'and_v3_id_fcomp.dat' using 1:2 with lines title 'x1,x2,x3,x4,x5', 'and_v3_id_fcomp.dat' using 1:3 with lines title 'x6', 'and_v3_id_fcomp.dat' using 1:4 with lines title 'x7'