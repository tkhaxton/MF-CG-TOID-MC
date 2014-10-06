
set logscale y
set format y "10^{%L}"

set key title ""

set xlabel "q_{xy} (A^{-1})"

filemonolayer="../../output/monolayer/post_analysis/xrd.allatom.600000-1000000.faceon.radav.total"
filebilayer="../../output/bilayer/post_analysis/xrd.allatom.600000-1000000.faceon.radav.total"
filestack="../../output/bilayer_stack_evaporate/post_analysis/xrd.allatom.500000-500000.faceon.radav.total"

plot filestack u 1:((column(2)+column(3)+column(4))) with lines lw 6 lt 1 title "stacked bilayer", \
filebilayer u 1:((column(2)+column(3))) with lines lw 6 lt 7 title "bilayer", \
filemonolayer u 1:((column(2))) with lines lw 6 lt 3 title "monolayer"

reset

