# *lsa*BGC-Pan

### *lsa*BGC-Pan - *mine the pan-BGC-ome of a microbial taxon for biosynthetic golden nuggets.*

*lsa*BGC-Pan reconfigures [*lsa*BGC](https://github.com/Kalan-Lab/lsaBGC) for easier installation via both Docker and Bioconda and features a new workflow bearing the same name as the repo. In addition to easier usability, there are some new analytical modules -e.g. (de-)association testing of BGC ortholog groups and GCFs and an improved framework for inferring horizontal transfer. Another key advantage is that it allows for joint analysis of both AntiSMASH and GECCO BGC predictions for a set of samples/genomes. 

The workflow has been tested to work with both bacterial and fungal genomes, but if you stumble on a bug/unexpected-error please open a GitHub issue and let us know!

<img src="https://github.com/Kalan-Lab/lsaBGC-Pan/assets/4260723/3aa3426e-39d6-4d25-91a3-44be288a6ad4" width="300">

### Manuscript:

> [Evolutionary investigations of the biosynthetic diversity in the skin microbiome using *lsa*BGC](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988). *Microbial Genomics 2023.* Rauf Salamzade, J.Z. Alex Cheong, Shelby Sandstrom, Mary Hannah Swaney, Reed M. Stubbendieck, Nicole Lane Starr, Cameron R. Currie, Anne Marie Singh, and Lindsay R. Kalan

### Documentation:

Documentation and three separate tutorial showing application to:

* A small set of 7 *Cutibacterium avidum* and *Cutibacterium acnes* genomes (test dataset included with repo)
* A set of 12 *Aspergillus flavus* genomes (fungal showcase)
* A set of 23 distinct *Streptomyces olivaceus* strains

can be found on the wiki at: https://github.com/Kalan-Lab/lsaBGC-Pan/wiki

### License:

```
BSD 3-Clause License

Copyright (c) 2024, Kalan-Lab

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
