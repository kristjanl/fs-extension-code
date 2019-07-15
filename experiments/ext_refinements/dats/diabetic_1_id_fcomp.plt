set term png
set output './dat_images/refs_diabetic_1.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:8]
plot 'diabetic_1_id_fcomp.dat' using 1:2 with lines title 'G', 'diabetic_1_id_fcomp.dat' using 1:3 with lines title 'X', 'diabetic_1_id_fcomp.dat' using 1:4 with lines title 'I', 'diabetic_1_id_fcomp.dat' using 1:5 with lines title 'T'