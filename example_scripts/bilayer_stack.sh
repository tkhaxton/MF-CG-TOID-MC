#!/bin/bash

inputfile="../output/bilayer/coord"
directory="../output/bilayer_stack"
mkdir -p $directory

time ../code/peptoid -inputfile $inputfile -base $directory -replicates 2 -ictype 7 -stackx 0 -stacky 0 -stackz 60 -cycles 100000 -wholebilayertranslatefreq 0.1



