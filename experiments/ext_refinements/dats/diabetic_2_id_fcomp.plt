set term png
set output './dat_images/refs_diabetic_2.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:6]
plot 'diabetic_2_id_fcomp.dat' using 1:2 with lines title 'G,X,I', 'diabetic_2_id_fcomp.dat' using 1:3 with lines title 'T'