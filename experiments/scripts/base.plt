set term png
$NAME$
set style line 1 linecolor rgb "blue"
set autoscale
unset label
set xtic auto
set ytic auto
set xlabel "t"
$VAR_LABEL$
plot '-' notitle with lines ls 1
$DATA$
e
