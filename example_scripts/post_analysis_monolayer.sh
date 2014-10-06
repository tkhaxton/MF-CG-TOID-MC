#!/bin/bash

trajectoryfile="../output/monolayer/trajectory"
directory="../output/monolayer/post_analysis"
mkdir -p $directory

time ../code/post_analysis -trajectoryfile $trajectoryfile -directory $directory -chainlength 28 -Nchains 48 -mincycle 500000 -nframes 5 -xrdcode 2



