# ScPathoQuant Installation and Usage

## Installation 
Install scPathoQuant repository via git clone:
```
git clone https://github.com/galelab/scPathoQuant.git
cd scPathoQuant
```

Using the condaenv.yml file, create conda environment. All required dependencies are included:

```
conda env create -f condaenv.yml -n scpathoquant
conda activate scpathoquant
pip install --no-deps .
```
ScPathoQuant can also be ran using Docker. For the purposes of this project, all analysis was done using conda when possible.


## Running scPathoQuant
ScPathoQuant required ```cellranger count``` output as input. Please refer to cellranger pre-processing documentation for more details. 

### Parameters
* -10x = path to cellranger count output
* -op = output directory
* -p = number of processors 
* -p2genome = path to reference genome sequence (and optionally pathogen gtf file)
* -align = alignment tool bbmap (default) or bowtie2
* --bbmap_params = bbmap specific parameters
* --bowtie2_params = bowtie2 specific parameters 

### General Example Run
```
scpathoquant -10x path/to/cellranger_output/ -op path/for/results/directory -p 8 -p2genome path/to/pathogen_reference_genome/ 
```
### Example Run: HIV Infected 1- 200k Reads Subset Aligned to Human Reference Genome
```
# run this in the ~/data directory
# run with default setting (bbmap as aligner)
scpathoquant -10x ./processed/SRR6825024-200k/ -op ./scpathoquant_outputs -p 8 -p2genome path/to/pathogen_reference_genome/
```

## Output - in progress
Outputs are in ~/data/scpathoquant_outputs folds 
