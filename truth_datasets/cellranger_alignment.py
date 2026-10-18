import argparse
import sys
import subprocess
from pathlib import Path
import time

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Running cellranger to align raw data to human genome with integrated pathogen genome. \n Adapted from analysis bash script from Whitmore et. al for use with Python.')
    parser.add_argument('-i', '--input',
        help='path to directory containing fastq files',
        required='True'
        )
    parser.add_argument('-g', '--genome_directory',
            help = 'path cell ranger reference genome',
            required = 'True')
    parser.add_argument('-o', '--output',
        help='output directory path',
        required='True'
        )
    return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output
genome_directory = args.genome_directory

# Cell Ranger resources
LOCAL_CORES = 5
LOCAL_MEM = 32

# parsing out how many samples will be processed
file_names = [file.name for file in Path(infile).iterdir() if file.is_file()]

# extracting the srrs
samples = set()
for file in file_names:
    sample_name = file.split('_')[0]
    samples.add(sample_name)

# making output directory and moving into it
subprocess.run(f'mkdir -p {outfile}', shell = True)
subprocess.run(f'cd {outfile}', shell = True)

print("========================================")
print("Cell Ranger processing")
print("========================================")
print(f"FASTQ directory: {infile}")
print(f'Reference: {genome_directory}')
print(f'Number of samples: {len(samples)}')
print(f'CPU cores: {LOCAL_CORES}')
print(f'Memory: {LOCAL_MEM} GB')
print('Processing one sample at a time')
print("========================================")

counter = 0
for sample_name in samples:
    counter += 1
    print("----------------------------------------")
    print(f"Processing sample: {sample_name}")
    print(f"Sample {counter} of {len(samples)}")
    print("----------------------------------------")

    # expected Cell Ranger FASTQ files 
    # hard coding the fastq read count for now, will find a better way to do this later -> only need below if I am doing error checking
    '''
    R1=f"{infile}/{sample_name}_200k_1.fastq"
    R2=f"{infile}/{sample_name}_200k_2.fastq"
    '''

    # running cell ranger
    cellranger_command = f'cellranger count --id {sample_name}--transcriptome {genome_directory} --fastqs {infile} --sample {sample_name} --localcores {LOCAL_CORES} --localmem {LOCAL_MEM} --create-bam false'

    start_time = time.perf_counter()
    result = subprocess.run(cellranger_command, shell = True)
    end_time = end_time = time.perf_counter()
    elapsed_time = start_time - end_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")

    status = result.returncode

    if status != 0:
        print()
        print(f"ERROR: Cell Ranger failed for {sample_name}")
        print(f"Exit status: {status}")
        sys.exit(status)

    print()
    print(f"Finished {sample_name}")
    print()

print('Finished processing all samples')