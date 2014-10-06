
filename="../../output/monolayer/timeseries"

set ylabel "Area/monomer (A^2)"
set xlabel "MC cycles"
cyclefactor=1
set logscale x
set format x "10^{%L}"
set xr [1e3:*]

plot filename u 1:(column(24)*column(25)/28/48) with lines lw 6 title ""

reset

