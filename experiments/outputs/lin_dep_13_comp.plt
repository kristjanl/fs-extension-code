set term png
set output './images/lin_dep_13_comp.png'
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x1"
plot '-' notitle with lines ls 1
e
