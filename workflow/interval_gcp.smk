import helperFun
import yaml
import pandas as pd
configfile: "config/config.yaml"
res_config = yaml.safe_load(open("config/resources.yaml"))

helperFun.make_temp_dir()
samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
org_ref = set(zip(samples.Organism.tolist(), samples.refGenome.tolist()))  # To get only unique combinations of organism and ref accession.
ORGANISM, REFGENOME = map(list, zip(*org_ref))  # Split above back to indiviudal lists for expand. There probably is a better way?

localrules: picard_intervals, create_intervals, index_ref
include: "rules/common.smk"

rule all:
    input:
        expand(config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "done.txt", Organism=ORGANISM, refGenome=REFGENOME),

rule index_ref:
    input:
        ref = ancient(config["refGenomeDir"] + "{refGenome}.fna")
    output: 
        indexes = expand(config["refGenomeDir"] + "{{refGenome}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"]),
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict"
    conda:
        "envs/gcp_mapping.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['index_ref']['mem']
    log:
        "logs/index_ref/{refGenome}.log" 
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} 2>> {log}
        samtools dict {input.ref} > {output.dictf} 2>> {log}
        """

rule picard_intervals:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
    output:
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list"    
    params:
        minNmer = int(config['minNmer'])
    conda:
        'envs/intervals.yml'
    log:
        "logs/{Organism}/{refGenome}/picard_intervals/log"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} &> {log}\n" 

rule create_intervals:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = ancient(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list")
    params:
        maxIntervalLen = int(config['maxIntervalLen']),
        maxBpPerList = int(config['maxBpPerList']),
        maxIntervalsPerList = int(config['maxIntervalsPerList']),
        minNmer = int(config['minNmer']),
        max_intervals = config['maxNumIntervals']
    output: 
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed",
        #config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "num_intervals.txt",
        l = directory(config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "lists/"),
        done = touch(config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "done.txt")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        Path(output.l).mkdir(exist_ok=True)
        if config['split_by_n']:
            NUM_OF_LISTS = helperFun.createListsGetIndices(params.maxIntervalLen, params.maxBpPerList, params.maxIntervalsPerList, params.minNmer, config["output"], config["intDir"], wildcards, input.dictf, input.intervals)
        else:
            NUM_OF_LISTS = make_intervals(output, config["output"], config["intDir"], wildcards, input.dictf, params.max_intervals)
