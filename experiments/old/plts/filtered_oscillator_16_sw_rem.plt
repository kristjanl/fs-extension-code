set term png
set output './images/filtered_oscillator_16_sw_rem_rem_QR_18_t_0.0.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x18"
set xrange [0.0:0.05]
set yrange [-2.6e-14:2.6e-14]
plot '-' notitle with lines ls 1
0 -0.000000000000008
0 0.000000000000008
0.05 0.000000000000008
0.05 -0.000000000000008
0 -0.000000000000008


e
