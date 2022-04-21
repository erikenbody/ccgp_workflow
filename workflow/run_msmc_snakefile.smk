import pandas as pd
import yaml
import helperFun
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

configfile: "config/config.yaml"
res_config = yaml.safe_load(open("config/resources.yaml"))

helperFun.make_temp_dir()
samples = pd.read_table(config["samples"], sep=",", dtype=str).replace(' ', '_', regex=True)
species_counts = samples.drop_duplicates(subset = ["BioSample", "refGenome", "Organism"]).value_counts(subset=['refGenome', 'Organism'])  #get BioSample for each refGenome/Organism combination
REFGENOME,ORGANISM = map(list, zip(*species_counts.index))  # split index into ref genome and organism

if config['remote_reads']:
    GS = GSRemoteProvider()
    GS_READS_PREFIX = config['remote_reads_prefix']

rule all:
    input:
        ancient(expand(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.pdf", zip, Organism=ORGANISM, refGenome=REFGENOME)),
        ancient(expand(config['output'] + "{Organism}/{refGenome}/" + "cleanup_confirmation.txt", zip, Organism=ORGANISM, refGenome=REFGENOME))

include: "rules/common.smk"
include: "rules/msmc.smk"