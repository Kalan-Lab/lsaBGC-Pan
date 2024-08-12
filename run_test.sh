#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz 
cd test_case/

# run lsaBGC on test set of select Cutibacterium avidium and Cutibacterium acnes genomes:
lsaBGC-Pan -g input_genomes/ -o lsaBGC-Pan_Results/ -c 4 -nb
