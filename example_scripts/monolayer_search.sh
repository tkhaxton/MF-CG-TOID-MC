#!/bin/bash

directory="../output/monolayer_search"
mkdir -p $directory

time ../code/monolayer -base $directory -chemistry 1 1 12 2 1 2 7 3 7 -interchainspacing 4.5 -temperature 10000 -cycles 10000



