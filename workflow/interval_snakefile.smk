import pandas as pd
import yaml
import helperFun
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

configfile: "config/config.yaml"
res_config = yaml.safe_load(open("config/resources.yaml"))

samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
species_counts = samples.drop_duplicates(subset = ["BioSample", "refGenome", "Organism"]).value_counts(subset=['refGenome', 'Organism'])  #get BioSample for each refGenome/Organism combination
REFGENOME,ORGANISM = map(list, zip(*species_counts.index))  # split index into ref genome and organism


rule all:
    input:
        ancient(expand(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed", zip, Organism=ORGANISM, refGenome=REFGENOME)),
        # ancient(expand(config['output'] + "{Organism}}/{refGenome}/" + config["intDir"] + "list{list}.list", zip, Organism=ORGANISM, refGenome=REFGENOME, list=config['maxNumIntervals'])),
        ancient(expand(config['output'] + "{Organism}/{refGenome}/" + "bam_sumstats.tsv", zip, Organism=ORGANISM, refGenome=REFGENOME)),
        expand(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "list{list}.list", Organism=ORGANISM, refGenome=REFGENOME, list=range(config['maxNumIntervals']))

include: "rules/common.smk"
include: "rules/intervals.smk"
include: "rules/msmc.smk"
include: "rules/sentieon.smk"
include: "rules/sumstats.smk"