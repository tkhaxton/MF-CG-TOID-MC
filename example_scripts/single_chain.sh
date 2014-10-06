#!/bin/bash

directory="../output/single_chain"
mkdir -p $directory

time ../code/peptoid -base $directory -chemistry 1 1 1 2 1 2 7 3 7 -interface 0 -ictype 1 -density 0.00001 -cycles 1000000



