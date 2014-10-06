

file="../../output/bilayer_stack_evaporate/timeseries"

set logscale x
set format x "10^{%L}"
set xr [1e3:*]

set xlabel "MC cycles"
set ylabel "Lamellar spacing (A)"

plot file u (column(1)):(column(34)/2) with lines lw 8 title ""

set term pop
set output
reset

