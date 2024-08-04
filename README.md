# *lsa*BGC-Pan

### *lsa*BGC-Pan - *mine the pan-BGC-ome of a microbial taxon for biosynthetic golden nuggets.*

*lsa*BGC-Pan reconfigures [*lsa*BGC](https://github.com/Kalan-Lab/lsaBGC) for easier installation via either Bioconda or Docker (_in progress_) and features a new workflow bearing the same name as the repo. In addition to easier usability, there are some new analytical modules -e.g. (de-)association testing of BGC ortholog groups and GCFs and an improved framework for inferring horizontal transfer. 

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

✔️ Tested Bioconda installation on both macOS & Linux with easy-to-use Docker image & wrapper script coming soon!

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
