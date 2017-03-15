set term png
set output './images/filtered_oscillator_32_plain.png
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x1"
plot '-' notitle with lines ls 1
0 0.2
0 0.35
0.05 0.35
0.05 0.2
0 0.2


e
