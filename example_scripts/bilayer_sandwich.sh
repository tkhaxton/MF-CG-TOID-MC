#!/bin/bash

inputfile="../output/monolayer_search/c10000.best.config"
directory="../output/bilayer_sandwich"
mkdir -p $directory

time ../code/bilayer_from_monolayer_enumerate -inputfile $inputfile -outputfile $directory/grid -chemistry 1 1 24 2 1 2 7 3 7 -firstx 0 -lastx 2 -xspace 0.2 -firsty 0 -lasty 2 -yspace 0.2 -firstz 7.5 -lastz 16 -zspace 0.5



