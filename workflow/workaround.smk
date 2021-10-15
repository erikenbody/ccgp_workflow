import os
from collections import defaultdict
import sys
import yaml
import helperFun
import pandas as pd
import glob
configfile: "config.yaml"
res_config = yaml.load(open("resources.yaml"))

bams = glob.glob("results/*/*/01_mappedReads/*.bam")
good_bams = [os.path.basename(bam).split("_final")[0] for bam in bams]
samples = pd.read_table(config["samples"], sep=",").replace(' ', '_', regex=True)
samples = samples[samples['BioSample'].isin(good_bams)]
org_ref = set(zip(samples.Organism.tolist(), samples.refGenome.tolist()))  # To get only unique combinations of organism and ref accession.
ORGANISM, REFGENOME = map(list, zip(*org_ref))  # Split above back to indiviudal lists for expand. There probably is a better way?
print(ORGANISM, REFGENOME)


rule all:
    input: 
        expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{lst}.raw.g.vcf.gz", Organism=ORGANISM, refGenome=REFGENOME, sample=samples.BioSample.tolist(), lst=list(range(1,31)))

rule bam2gvcf:
    """
    This rule scatters analyses over two dimensions: sample name and list file. For each BAM file, one per sample,
    a GVCF is created for all the scaffolds present in a given list file.
    """
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        bam = config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "{sample}" + config['bam_suffix'],
        l = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "list0.list",
        gvcf = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L0.raw.g.vcf.gz",
        gvcf_idx = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L0.raw.g.vcf.gz.tbi",
    output: 
        gvcf = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz",
        gvcf_idx = config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz.tbi",
        doneFile = touch(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.done")
    resources: 
        #!The -Xmx value the tool is run with should be less than the total amount of physical memory available by at least a few GB
        # subtract that memory here
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2gvcf']['mem'],   # this is the overall memory requested
        reduced = lambda wildcards, attempt: attempt * (res_config['bam2gvcf']['mem'] - 3000)  # this is the maximum amount given to java

    params:
        minPrun = config['minP'],
        minDang = config['minD']
    shell:
        """
        cp {input.gvcf} {output.gvcf}
        cp {input.gvcf_idx} {output.gvcf_idx}
        """




