set term png
set output './dat_images/refs_lin_dep_20.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [4:7]
plot 'lin_dep_20_id_fcomp.dat' using 1:2 with lines title 'x1', 'lin_dep_20_id_fcomp.dat' using 1:3 with lines title 'x2', 'lin_dep_20_id_fcomp.dat' using 1:4 with lines title 'x3', 'lin_dep_20_id_fcomp.dat' using 1:5 with lines title 'x4', 'lin_dep_20_id_fcomp.dat' using 1:6 with lines title 'x5', 'lin_dep_20_id_fcomp.dat' using 1:7 with lines title 'x6', 'lin_dep_20_id_fcomp.dat' using 1:8 with lines title 'x7', 'lin_dep_20_id_fcomp.dat' using 1:9 with lines title 'x8', 'lin_dep_20_id_fcomp.dat' using 1:10 with lines title 'x9', 'lin_dep_20_id_fcomp.dat' using 1:11 with lines title 'x10', 'lin_dep_20_id_fcomp.dat' using 1:12 with lines title 'x11', 'lin_dep_20_id_fcomp.dat' using 1:13 with lines title 'x12', 'lin_dep_20_id_fcomp.dat' using 1:14 with lines title 'x13', 'lin_dep_20_id_fcomp.dat' using 1:15 with lines title 'x14', 'lin_dep_20_id_fcomp.dat' using 1:16 with lines title 'x15', 'lin_dep_20_id_fcomp.dat' using 1:17 with lines title 'x16', 'lin_dep_20_id_fcomp.dat' using 1:18 with lines title 'x17', 'lin_dep_20_id_fcomp.dat' using 1:19 with lines title 'x18', 'lin_dep_20_id_fcomp.dat' using 1:20 with lines title 'x19', 'lin_dep_20_id_fcomp.dat' using 1:21 with lines title 'x20'