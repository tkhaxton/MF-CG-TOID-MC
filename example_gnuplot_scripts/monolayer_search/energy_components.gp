
filename="../../output/monolayer_search/c10000.best.energy"

backbonebulk=-10.1
phenylbulk=-0.76
aminobulk=-71.3
carboxylbulk=-79.9
phenylfrac=0.5
aminofrac=0.25
carboxylfrac=0.25
offset=backbonebulk+phenylbulk*phenylfrac+aminobulk*aminofrac+carboxylbulk*carboxylfrac

set ylabel "E/monomer (kcal/mol)"
set xlabel "MC cycles"
set logscale x
set format x "10^{%L}"
set xr [*:*]

set key left Left reverse bottom

plot filename u 1:(column(12)+column(13)+column(14)+column(15)+column(16)+column(17)+column(18)+column(19)-offset) with lines lw 6 title "solvation", \
filename u 1:2 with lines lw 6 title "backbone", \
filename u 1:3 with lines lw 6 title "sidechain", \
filename u 1:(column(4)+column(8)) with lines lw 6 title "nonpolar", \
filename u 1:(column(6)+column(7)+column(10)+column(11)) with lines lw 6 title "charged", \
filename u 1:(column(2)+column(3)+column(4)+column(5)+column(6)+column(7)+column(8)+column(9)+column(10)+column(11)+column(12)+column(13)+column(14)+column(15)+column(16)+column(17)+column(18)+column(19)-offset) with lines lw 6 title "total"
reset


