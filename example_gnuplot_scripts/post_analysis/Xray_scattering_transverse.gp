
set logscale y
set format y "10^{%L}"

set key title ""

set xlabel "q_{z} (A^{-1})"

filestack="../../output/bilayer_stack_evaporate/post_analysis/xrd.allatom.500000-500000.transverse.total"

plot filestack u 1:((column(2)+column(3)+column(4))) with lines lw 6 lt 1 title "stacked bilayer"

reset

