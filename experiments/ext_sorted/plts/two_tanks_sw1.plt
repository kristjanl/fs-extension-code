set term png
set output './images/two_tanks_sw1_method_1_t_1.99.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x1"
set xrange [0.0:2.2]
set yrange [-57.2913810353:56.0522082583]
plot '-' notitle with lines ls 1
0 1.5
0 2.5
0 2.5
0 1.5
0 1.5


0 1.5
0 2.5
0.01 2.5
0.01 1.5
0 1.5


e
