

file="../../output/bilayer/timeseries"

set logscale x
set format x "10^{%L}"
set xr [1e3:*]

set xlabel "MC cycles"
set ylabel "Area/monomer (A^2)"

plot file u (column(1)):(column(28)*column(29)/28/48) with lines lw 8 title ""

set term pop
set output
reset

