set term png
set output './images/bouncing_ball_sw_10.png'
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x1"
plot '-' notitle with lines ls 1
0.000000000001 9.95095
0.1 9.95095
0.1 10.2
0.000000000001 10.2
0.000000000001 9.95095


e
