#!/bin/bash

inputfile="../output/bilayer_search/c10000.best.config"
directory="../output/bilayer"
mkdir -p $directory

time ../code/peptoid -inputfile $inputfile -base $directory -chemistry 1 1 96 2 1 2 7 3 7 -xchains 2 -interface 0 -ictype 3 -cycles 1000000



