# ScPathoQuant Installation and Usage

## Installation 

Install scPathoQuant repository via git clone and add to your path:
```
git clone https://github.com/galelab/scPathoQuant.git
cd scPathoQuant
```
Build ScPathoQuant environment
```
conda env create -f envs/scpathoquant.yml
```


## Running scPathoQuant
ScPathoQuant requires ```cellranger count``` output as input. Please refer to cellranger pre-processing documentation for more details. 

### Parameters
For the purpose of this analysis, bbmap was the only alignment only alignment tool investigated because it is splice-aware.

* -10x = path to cellranger count output
* -op = output directory
* -p = number of processors 
* -p2genome = path to reference genome sequence (and optionally pathogen gtf file)
* -align = alignment tool bbmap (default) or bowtie2
* --bbmap_params = bbmap specific parameters 
* --bowtie2_params = bowtie2 specific parameters 

### General Example Run
```
scpathoquant -10x path/to/cellranger_output/ -op path/to/output/directory -p 8 -p2genome path/to/pathogen_reference_genome/ --bbmap_params="--semiperfectmode"
```

### Example Run: HIV Infected 1- 25M Reads Subset Aligned to Human Reference Genome
```
# processed cellranger count output can be accessed in the /home/skyemegan directory
scpathoquant -10x /home/skyemegan/SRR6825024-25M -op ./scpathoquant_outputs -p 8 -p2genome data/hiv_genome --bbmap_params="--semiperfectmode"
```

## Output 
For the purposes of this analysis, the only output of interest are:
* pathogen_al_counts.csv:  total number of reads mapping to the pathogen in each cell
