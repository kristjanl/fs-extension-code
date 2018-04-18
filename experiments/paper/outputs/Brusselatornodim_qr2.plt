set term png
set output './images/Brusselatornodim_qr2.png'
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
set ylabel "x1"
plot '-' notitle with lines ls 1
e
