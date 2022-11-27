# RNAseq

Trimming to featurecounts snakemake

Currently Loaded Modules:

`module load gbc-samtools/1.12 
gbc-hisat2/2.2.1 
gbc-cutadapt/1.16 
python/py37-anaconda-2019.10 
snakemake/5.7.1-py37 
gbc-fastqc 
gbc-subread`

Step-by-step of install and analysis

1. git clone this repository

`git clone  https://github.com/achisha21/RNAseq-snakemake.git`

2. Activate the python anaconda environment

`conda activate snakemake`

3. Edit the config.json file and cluster.json files

4. Ensure meta-data table contains all of the necessairy fields

** NOTE EXACT HEADERS HAVE TO BE ENFORCED or key errors will be thrown during processing**

5. Launch jobs
The use of --latency-wait allows for SLURM to catch up writing the files and posting the file handles so Snakemake can see them.

`snakemake --latency-wait 120 -p -j 100 --profile slurm`


