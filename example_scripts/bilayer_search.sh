#!/bin/bash

inputfile="../output/monolayer_search/c10000.best.config"
directory="../output/bilayer_search"
mkdir -p $directory

time ../code/bilayer_from_monolayer -inputfile $inputfile -base $directory -chemistry 1 1 24 2 1 2 7 3 7 -leafxoffsetfrac 0 -leafyoffsetfrac 0.6 -leafspacing 14.5 -temperature 10000 -cycles 10000



