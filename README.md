# *lsa*BGC-Pan
[![Documentation](https://img.shields.io/badge/Documentation-Wiki-darkgreen?style=flat-square&maxAge=2678400)](https://github.com/Kalan-Lab/lsaBGC-Pan/wiki)
[![Documentation](https://img.shields.io/badge/Installation-limegreen?style=flat-square&maxAge=2678400)](https://github.com/Kalan-Lab/lsaBGC-Pan/wiki)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lsabgc/badges/version.svg)](https://anaconda.org/bioconda/lsabgc)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lsabgc/badges/platforms.svg)](https://anaconda.org/bioconda/lsabgc)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lsabgc/badges/latest_release_date.svg)](https://anaconda.org/bioconda/lsabgc)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lsabgc/badges/downloads.svg)](https://anaconda.org/bioconda/lsabgc)
[![Docker](https://img.shields.io/badge/Docker-DockerHub-darkred?style=flat-square&maxAge=2678400)](https://hub.docker.com/r/raufs/lsabgc_pan)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/lsabgc/badges/license.svg)](https://anaconda.org/bioconda/lsabgc)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.13309206.svg)](https://doi.org/10.5281/zenodo.13309206)
[![Manuscript](https://img.shields.io/badge/Manuscript-MGen-darkblue?style=flat-square&maxAge=2678400)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988)

### *lsa*BGC-Pan - *mine the pan-BGC-ome of a microbial taxon for biosynthetic golden nuggets.*

*lsa*BGC-Pan reconfigures [*lsa*BGC](https://github.com/Kalan-Lab/lsaBGC) for easier installation via either Bioconda or Docker and features a new workflow bearing the same name as the repo. In addition to easier usability, there are some new analytical modules -e.g. (de-)association testing of BGC ortholog groups and GCFs and an improved framework for inferring horizontal transfer. 

<p align="center">
<img src="https://github.com/Kalan-Lab/lsaBGC-Pan/assets/4260723/3aa3426e-39d6-4d25-91a3-44be288a6ad4" width="300">
</p>

### Manuscript:

> [Evolutionary investigations of the biosynthetic diversity in the skin microbiome using *lsa*BGC](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988). *Microbial Genomics 2023.* Rauf Salamzade, J.Z. Alex Cheong, Shelby Sandstrom, Mary Hannah Swaney, Reed M. Stubbendieck, Nicole Lane Starr, Cameron R. Currie, Anne Marie Singh, and Lindsay R. Kalan

### Key Highlights:

✔️ Works for both fungi & bacteria

✔️ Allows for joint analysis of antiSMASH & GECCO BGC predictions (in analysis of _Streptomyces olivaceus_ - this leads to a 42.4% increase in distinct GCFs to using antiSMASH alone)

✔️ Better consideration for incomplete BGCs due to assembly fragmentation

✔️ New analytical features including: (1) genome-wide association testing of orthogroups with GCF co-occurence & (2) improved assessment of horizontal transfer for BGC-associated orthogroups

✔️ Improved consolidated spreadsheet that is easier to assess

✔️ Support for small scale (< 30 genomes) analysis on laptop with minimal databases (~5 GB).

✔️ Bioconda installation tested on both macOS & Linux with easy-to-use Docker image & wrapper script also coming soon!

✔️ BSD-3 License & no uploading of data to webservers = support for industry research

### Documentation:

Documentation and three separate tutorials showing application to:

* A small set of 7 *Cutibacterium avidum* and *Cutibacterium acnes* genomes (test dataset included with repo)
* A set of 12 *Aspergillus flavus* genomes (fungal showcase)
* A set of 23 distinct *Streptomyces olivaceus* strains (shows joint AntiSMASH + GECCO analysis)

can be found on the wiki at: https://github.com/Kalan-Lab/lsaBGC-Pan/wiki

### Example Commands:

Perform analysis using a directory of AntiSMASH results as input:

```bash
lsaBGC-Pan -a AntiSMASH_Results/ -o Pan_Results/ -c 10
```

Provide a directory of AntiSMASH results as input and incorporate GECCO BGC predictions as well:

```bash
lsaBGC-Pan -a AntiSMASH_Results/ -o Pan_Results/ -c 10 -rg
```

Provide a directory of genomes in FASTA format for GECCO-based BGC predictions and analysis (only works for bacteria): 

```bash
lsaBGC-Pan -g Directory_of_Genomes_in_FASTA/ -o Pan_Results/ -c 10
```
