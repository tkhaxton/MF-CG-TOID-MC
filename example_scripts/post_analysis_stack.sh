#!/bin/bash

trajectoryfile="../output/bilayer_stack_evaporate/trajectory"
directory="../output/bilayer_stack_evaporate/post_analysis"
mkdir -p $directory

time ../code/post_analysis -trajectoryfile $trajectoryfile -directory $directory -chainlength 28 -Nchains 192 -mincycle 400000 -nframes 1 -xrdcode 3



