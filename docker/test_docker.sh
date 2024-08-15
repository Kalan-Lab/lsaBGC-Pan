#!/usr/bin/env bash

#############################################################################################################
##### This is same test set as for conda based installation but using the run_lsaBGC-Pan.sh wrapper instead.
#############################################################################################################

# RUN THIS IN THE SAME DIRECTORY AS THE run_lsaBGC-Pan.sh SCRIPT

# change permissions for run_lsaBGC-Pan.sh in this directory to allow execution.
chmod +x run_lsaBGC-Pan.sh

# Step 0: Download and uncompress test_case.tar.gz and cd into it.
rm test_case.tar.gz
wget https://github.com/Kalan-Lab/lsaBGC-Pan/raw/main/test_case.tar.gz
rm -rf test_case/
tar -zxvf test_case.tar.gz
cd test_case/

ln -s ../run_lsaBGC-Pan.sh .

# run test
./run_lsaBGC-Pan.sh -i input_genomes/ -o lsabgc_pan_results/ -nb -c 4
