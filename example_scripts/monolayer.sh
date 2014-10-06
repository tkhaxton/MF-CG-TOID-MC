#!/bin/bash

inputfile="../output/monolayer_search/c10000.best.config"
directory="../output/monolayer"
mkdir -p $directory

time ../code/peptoid -inputfile $inputfile -base $directory -chemistry 1 1 48 2 1 2 7 3 7 -xchains 2 -interface 1 -ictype 2 -surfacepressure 0.044 -cycles 1000000



