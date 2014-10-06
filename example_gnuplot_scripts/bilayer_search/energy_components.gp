
filename="../../output/bilayer_search/c10000.best.energy"

set ylabel "E/monomer (kcal/mol)"
set xlabel "MC cycles"
set logscale x
set format x "10^{%L}"

set key left Left reverse bottom

plot filename u 1:2 with lines lw 6 title "backbone", \
filename u 1:3 with lines lw 6 title "sidechain", \
filename u 1:(column(4)+column(8)) with lines lw 6 title "nonpolar", \
filename u 1:(column(6)+column(7)+column(10)+column(11)) with lines lw 6 title "charged", \
filename u 1:(column(2)+column(3)+column(4)+column(5)+column(6)+column(7)+column(8)+column(9)+column(10)+column(11)) with lines lw 6 title "total"
reset


