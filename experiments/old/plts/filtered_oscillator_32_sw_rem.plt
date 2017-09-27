set term png
set output './images/filtered_oscillator_32_sw_rem_plain_rem_34_t_0.0.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x34"
set xrange [0.0:0.05]
set yrange [-2e-15:2e-15]
plot '-' notitle with lines ls 1
0 -0.000000000000002
0 0.000000000000002
0.05 0.000000000000002
0.05 -0.000000000000002
0 -0.000000000000002


e
