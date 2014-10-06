#!/bin/bash

inputfile="../output/bilayer_stack/coord"
directory="../output/bilayer_stack_evaporate"
mkdir -p $directory

time ../code/peptoid -inputfile $inputfile -base $directory -ictype 0 -cycles 500000 -wholebilayertranslatefreq 0.1 -shiftbilayergapfreq 0.1 -reset 1



