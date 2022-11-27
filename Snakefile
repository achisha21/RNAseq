#! /usr/bin/py
import pandas as pd 
import os

#trim -> fastqc -> mapping -> sort bam files -> index bam files -> featurecounts

configfile: "config.json"
localrules: all, mkdir 

df = pd.read_csv(config["meta_file"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)
print(df.index)

def get_pair_gz(sample_id):
    dir=config["raw_fastq_gz_dir"]
    return tuple(os.path.join(dir, df.loc[str(sample_id), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

def get_forward_primer(sample_id):
    return df.loc[sample_id]["Adapter_1"]

def get_reverse_primer(sample_id):
    return df.loc[sample_id]["Adapter_2"]

rule all:
    input:expand("{dir}/readCounts.tsv", dir=config["dir_names"]["readCounts_dir"], sample_id=sample_ids)
    run:
        for sample in sample_ids:
            print("Wrapping up pipeline")

rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs = list(config["dir_names"].values())
    resources: time_min=10, mem_mb=2000, cpus=1
    shell: "mkdir -p {params.dirs}"

rule trim: 
    input:
        rules.mkdir.output,
        all_read1 = lambda wildcards: get_pair_gz(wildcards.sample_id)[0],
        all_read2 = lambda wildcards: get_pair_gz(wildcards.sample_id)[1]
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        trimmed_read1 = config["dir_names"]["trimmed_dir"]+ "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_read2 = config["dir_names"]["trimmed_dir"]+ "/{sample_id}.trimmed.R2.fastq.gz",
        trimmed_stats = config["dir_names"]["trimmed_dir"]+ "/{sample_id}.trimmed.stats"
    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wildcards: get_forward_primer(wildcards.sample_id),
        adapter2= lambda wildcards: get_reverse_primer(wildcards.sample_id)
    shell: "cutadapt -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} >{output.trimmed_stats}"

rule fastqc:
    input:
        read1= rules.trim.output.trimmed_read1,
        read2= rules.trim.output.trimmed_read2
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        fastqc_file1 = config["dir_names"]["fastqc_dir"]+ "/{sample_id}.trimmed.R1_fastqc.html",
        fastqc_file2 = config["dir_names"]["fastqc_dir"]+ "/{sample_id}.trimmed.R2_fastqc.html"
    shell: "fastqc {input.read1} -o outputs/fastqc ; fastqc {input.read2} -q -o outputs/fastqc/" 

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.trim.output.trimmed_read2,
        p3 = rules.fastqc.output.fastqc_file1,
        p4 = rules.fastqc.output.fastqc_file2
    resources: time_min=360, mem_mb=20000, cpus=6
    output: 
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
        stats = config["dir_names"]["mapped_dir"]+"/{sample_id}.stats"
    params:
        thread = config["params"]["hisat2"]["threads"],
        map_all = config["params"]["hisat2"]["all"],
        reference = config["params"]["hisat2"]["hisat2_reference"]
    shell:
        """
        hisat2 --threads {params.thread} -1 {input.p1} -2 {input.p2} -x {params.reference} 2> {output.stats} | samtools sort -T {output.mapped_bam_file}.tmp -O bam -o {output.mapped_bam_file}
        """

rule sort_bam:
    input:
        sorted_bam = rules.map.output.mapped_bam_file
    output:
        sorted_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.sorted.bam"
    resources: time_min=360, mem_mb=20000, cpus=6
    shell: "samtools view -hb {input.sorted_bam} | samtools sort -T {input.sorted_bam}.tmp -o {output.sorted_bam_file}" 

rule index_bam:
    input:
        index_bam = rules.sort_bam.output.sorted_bam_file
    resources: time_min=360, mem_mb=2000, cpus=1
    output:
        index_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam.bai"
    shell: "samtools index {input.index_bam} {output.index_bam_file}"

rule featurecount:
    input:
        sorted_bams = rules.sort_bam.output.sorted_bam_file,
        index_bams = rules.index_bam.output.index_bam_file
    resources: time_min=360, mem_mb=20000, cpus=6
    output:
        counts_file = config["dir_names"]["readCounts_dir"] + "/readCounts.tsv" 
    shell: "featureCounts -a /projects/academic/gbcstaff/igenome-041714/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o {output.counts_file} -T 6 /projects/academic/gbcstaff/Intern_Projects/achisha/Blumental-HHYTHDRX2/outputs/mapping/*.sorted.bam"

