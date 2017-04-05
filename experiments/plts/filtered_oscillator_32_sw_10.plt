set term png
set output './images/filtered_oscillator_32_sw_10_9_t_0.0.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x9"
set xrange [0.0:0.05]
set yrange [-1e-10:3.7e-09]
plot '-' notitle with lines ls 1
0 -0.0000000001
0 0.0000000037
0.05 0.0000000037
0.05 -0.0000000001
0 -0.0000000001


e
