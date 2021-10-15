import glob
from os import path
localrules: collect_sumstats, download_reference

ruleorder: index_ref > download_reference
### RULES ###

checkpoint get_fastq_pe:
    output:
        temp(directory(config["fastqDir"] + "{Organism}/{sample}/{run}"))
        
    params:
        outdir = config["fastqDir"] + "{Organism}/{sample}/{run}",
        tmpdir = config['tmp_dir']
    conda: "../envs/fastq2bam.yml"
    threads: int(res_config['get_fastq_pe']['threads'])
    benchmark:
        "benchmarks/{Organism}/fasterq_dump/{sample}/{run}.log"
    log:
        "logs/{Organism}/fasterq_dump/{sample}/{run}.log"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['get_fastq_pe']['mem']
    shell:
        # need to change exit code to catch if exit code is not 0. should just use ||. then use keep going snakemake option
        """
        fasterq-dump {wildcards.run} -O {params.outdir} -t {params.tmpdir} -e {threads} &> {log}
        """


rule download_reference:
    output:
        outdir = directory(config["refGenomeDir"] + "{refGenome}"),
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        dataset = config["refGenomeDir"] + "{refGenome}_dataset.zip"
    log:
        "logs/dl_reference/{refGenome}.log"
    conda:
        "../envs/fastq2bam.yml"
    shell:
        "datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {output.dataset} {wildcards.refGenome} > {log}"
        "&& 7z x {output.dataset} -aoa -o{output.outdir}"
        "&& cat {output.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref}"


rule index_ref:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        dataset = config["refGenomeDir"] + "{refGenome}_dataset.zip"
    output: 
        indexes = expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"])
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['index_ref']['mem']
    log:
        "logs/index_ref/{refGenome}.log" 
    shell:
        "bwa index {input.ref} 2> {log}"

rule fastp_pe:
    input:
        config["fastqDir"] + "{Organism}/{sample}/{run}"
    output: 
        r1 = temp(config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}/{run}_1.fastq.gz"),
        r2 = temp(config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}/{run}_2.fastq.gz"),
    params:
        r1 = config["fastqDir"] + "{Organism}/{sample}/{run}/" + "{run}_1.fastq",
        r2 = config["fastqDir"] + "{Organism}/{sample}/{run}/" + "{run}_2.fastq"
    conda:
        "../envs/fastq2bam.yml"
    threads: res_config['fastp']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['fastp']['mem'] 
    shell:
        "fastp --in1 {params.r1} --in2 {params.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--thread {threads} "
        "--detect_adapter_for_pe"
rule fastp_se:
    input:
        config["fastqDir"] + "{Organism}/{sample}/{run}"
    output: 
        r1 = temp(config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}/{run}.fastq.gz"),
    params:
        r1 = config["fastqDir"] + "{Organism}/{sample}/{run}/" + "{run}.fastq",
    conda:
        "../envs/fastq2bam.yml"
    threads: res_config['fastp']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['fastp']['mem'] 
    shell:
        "fastp --in1 {params.r1} "
        "--out1 {output.r1} "
        "--thread {threads} "

rule bwa_map:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        reads = get_reads,
        # the following files are bwa index files that aren't directly input into command below, but needed
        indexes = expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"])
    output: 
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "preMerge/{sample}/{run}.bam")
    params:
        get_read_group
    benchmark:
        "benchmarks/{Organism}/{refGenome}/bwamap/{sample}/{run}.log"
    conda:
        "../envs/fastq2bam.yml"
    threads: res_config['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bwa_map']['mem'] 
    log:
        "logs/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    benchmark:
        "benchmarks/{Organism}/bwa/{refGenome}_{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {threads} {params} {input.ref} {input.reads} 2> {log} | samtools sort -o {output.bam} -"

rule merge_bams:
    input: lambda wildcards: expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    output: 
        bam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam"),
        bai = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam.bai")

    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['merge_bams']['mem']
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam}"

rule dedup:
    input: 
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam",
        bai = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam.bai"
    output:
        dedupBam = temp(config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix']),
        dedupMet = temp(config['output'] + "{Organism}/{refGenome}/" + config['sumstatDir'] + "{sample}_dedupMetrics.txt")
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['dedup']['mem'] 
    log:
        "logs/{Organism}/dedup/{refGenome}_{sample}.txt"
    benchmark:
        "benchmarks/{Organism}/dedup/{refGenome}_{sample}.txt"
    shell:
        "picard MarkDuplicates I={input[0]} O={output.dedupBam} METRICS_FILE={output.dedupMet} REMOVE_DUPLICATES=false TAGGING_POLICY=All &> {log}\n"
        "picard BuildBamIndex I={output.dedupBam} &> {log}"
