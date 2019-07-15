set term png
set output './dat_images/refs_bouncing_ball.png'
set datafile separator ';'
set xlabel 'step'
set ylabel 'iters'
set autoscale
set yrange [0:2]
plot 'bouncing_ball_id_fcomp.dat' using 1:2 with lines title 'x', 'bouncing_ball_id_fcomp.dat' using 1:3 with lines title 'v'