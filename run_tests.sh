#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz 
cd test_case/

# Step 1: create genome listing inputs for lsaBGC-Ready.py
lsaBGC-Pan -i Primary_Genomes/ 

lsaBGC-Pan -a Primary_Genome_AntiSMASH_Results/ 