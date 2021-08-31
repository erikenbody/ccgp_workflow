import os
from collections import defaultdict
import sys
import yaml
import helperFun
import pandas as pd
import glob
configfile: "config/config.yaml"
res_config = yaml.load(open("config/resources.yaml"))

bams = glob.glob("results/*/*/01_mappedReads/*.bam")
good_bams = [os.path.basename(bam).split("_final")[0] for bam in bams]
samples = pd.read_table(config["samples"], sep=",").replace(' ', '_', regex=True)
samples = samples[samples['BioSample'].isin(good_bams)]
org_ref = set(zip(samples.Organism.tolist(), samples.refGenome.tolist()))  # To get only unique combinations of organism and ref accession.
ORGANISM, REFGENOME = map(list, zip(*org_ref))  # Split above back to indiviudal lists for expand. There probably is a better way?
print(ORGANISM, REFGENOME)
### workflow ###

rule all:
    input:
        expand(config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}.final.vcf.gz", zip, Organism=ORGANISM, refGenome=REFGENOME)

include: "rules/common.smk"
include: "rules/bam2vcf_gatk.smk"




