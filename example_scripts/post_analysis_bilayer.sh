#!/bin/bash

trajectoryfile="../output/bilayer/trajectory"
directory="../output/bilayer/post_analysis"
mkdir -p $directory

time ../code/post_analysis -trajectoryfile $trajectoryfile -directory $directory -chainlength 28 -Nchains 96 -mincycle 500000 -nframes 5 -xrdcode 2



