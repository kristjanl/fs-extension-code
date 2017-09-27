set term png
set output './images/rod_reactor_sw_10step_infl_3_t_0.9.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x3"
set xrange [0.0:1.0]
set yrange [19.99999:21.00001]
plot '-' notitle with lines ls 1
0 20
0 20.1
0.1 20.1
0.1 20
0 20


0.1 20.1
0.1 20.2
0.2 20.2
0.2 20.1
0.1 20.1


0.2 20.2
0.2 20.3
0.3 20.3
0.3 20.2
0.2 20.2


0.3 20.3
0.3 20.4
0.4 20.4
0.4 20.3
0.3 20.3


0.4 20.4
0.4 20.5
0.5 20.5
0.5 20.4
0.4 20.4


0.5 20.5
0.5 20.6
0.6 20.6
0.6 20.5
0.5 20.5


0.6 20.6
0.6 20.7
0.7 20.7
0.7 20.6
0.6 20.6


0.7 20.7
0.7 20.8
0.8 20.8
0.8 20.7
0.7 20.7


0.8 20.8
0.8 20.9
0.9 20.9
0.9 20.8
0.8 20.8


0.9 20.9
0.9 21
1 21
1 20.9
0.9 20.9


e
