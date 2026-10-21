# Viral-Track Installation and Usage
**Current Progress**: Viral-Track is installed but not running, as pre-processing steos are still in progress.

## Dependencies
* Spliced Transcript Alignment to A Reference (STAR) (http://github.com/alexdobin/STAR)
* Samtools (https://www.htslib.org/download/)
* StringTie (https://ccb.jhu.edu/software/stringtie/)
* UMI-tools (https://github.com/CGATOxford/UMI-tools)

**Next Step**: creating a conda environment for these dependencies for ease of access

## Installation 
Viral-Track is a collection of R scripts rather than a stand-alone tool. These separate scripts can be accessed via git clone.
```
git clone https://github.com/PierreBSC/Viral-Track.git
cd Viral-Track
```
**Next step**: create one script/command that will call each script as needed (vs running multiple different scripts)

## Data Pre-processing
### Index and Annotation File
Create STAR index for both host and virus reference genomes
