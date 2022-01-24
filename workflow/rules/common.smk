import glob
import re
import os
from pathlib import Path
from collections import defaultdict, deque
from snakemake.exceptions import WorkflowError
### INPUT FUNCTIONS ###

def get_gvcfs(wildcards):
    sample_names = samples.BioSample.tolist()
    return expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz", **wildcards, sample=sample_names)


def get_bams_for_dedup(wildcards):
    runs = samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist()
    
    if len(runs) == 1:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=runs)
    else:
        return config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam"
def get_bai_for_dedup(wildcards):
    runs = samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist()
    if len(runs) == 1:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam.bai", run=runs)
    else:
        return config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam.bai"
def get_bams_for_merge(wildcards):
    runs = samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist()
    if len(runs) == 1:
        return "non-existing-filename"  # this is a hack to coerce snakemake to not execute the rule this feeds if there is only one bam file
    else:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=runs)

def get_reads(wildcards):
    """Returns local read files if present. Defaults to SRR if no local reads in sample sheet."""
    row = samples.loc[samples['Run'] == wildcards.run]
    if 'fq1' in samples.columns and 'fq2' in samples.columns:
        #if os.path.exists(row.fq1.item()) and os.path.exists(row.fq2.item()):
        r1 = row.fq1.item()
        r2 = row.fq2.item()
        return {"r1": r1, "r2": r2}
    #    else:
    #        print(row.fq1.item(), row.fq2.item())
    #        raise WorkflowError(f"fq1 and fq2 specified for {wildcards.sample}, but files were not found.")
    else:
        r1 = config["fastqDir"] + f"{wildcards.Organism}/{wildcards.sample}/{wildcards.run}_1.fastq.gz",
        r2 = config["fastqDir"] + f"{wildcards.Organism}/{wildcards.sample}/{wildcards.run}_2.fastq.gz"
        return {"r1": r1, "r2": r2}
def get_read_group(wildcards):
    """Denote sample name and library_id in read group."""
    return r"-R '@RG\tID:{lib}\tSM:{sample}\tPL:ILLUMINA'".format(
        sample=wildcards.sample,
        lib=samples.loc[samples['BioSample'] == wildcards.sample]["LibraryName"].tolist()[0]
    )

def get_sumstats(wildcards):
    # Gets the correct sample given the organism and reference genome
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    fastpFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_fastp.out", sample=_samples)
    dedupFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_dedupMetrics.txt", sample=_samples)
    alnSumMetsFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_AlnSumMets.txt", sample=_samples)
    coverageFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_coverage.txt", sample=_samples)
    validateFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_validate.txt", sample=_samples)
    return {'fastpFiles': fastpFiles, 'dedupFiles': dedupFiles, 'alnSumMetsFiles': alnSumMetsFiles, 'coverageFiles': coverageFiles, 'validateFiles': validateFiles}

def get_gather_vcfs(wildcards):
    """
    Gets filtered vcfs for gathering step. This function gets the interval list indicies from the corresponding
    genome, then produces the file names for the filtered vcf with list index."""
    
    #checkpoint_output = checkpoints.create_intervals.get(**wildcards).output[0]
    
    
    path = Path(workflow.default_remote_prefix, config['output'], wildcards.Organism,wildcards.refGenome, config['intDir'], "lists/")
    print(path)
    lists = range(len(list(path.glob("*.list"))))
    return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['vcfDir_gatk'] + "filtered_L{index}.vcf", index=lists)

def gather_vcfs_CLI(wildcards):
    """
    Gatk enforces that you have a -I before each input vcf, so this function makes that string
    """
    vcfs = expand(config['output'] + "{Organism}/{refGenome}/" + config['vcfDir_gatk'] + "filtered_L{index}.vcf", **wildcards, index=range(10))
    #print(vcfs)
    
    out = " ".join(["-I " + vcf for vcf in vcfs])
    return out

def write_db_mapfile(wildcards, output):
    dbMapFile = output.dbMapFile
    sample_names = set(samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist())
    with open(dbMapFile, 'w') as f:
        for sample in sample_names:
            gvcf_path = os.path.join(workflow.default_remote_prefix,config['output'], wildcards.Organism, wildcards.refGenome, config['gvcfDir'], sample, f"L{wildcards.list}.raw.g.vcf.gz") 
            print(sample, gvcf_path, sep="\t", file=f)

def get_input_for_mapfile(wildcards):
    sample_names = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    gvcfs = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz", sample=sample_names, **wildcards)
    gvcfs_idx = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz.tbi", sample=sample_names, **wildcards)
    doneFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.done", sample=sample_names, **wildcards)
    return gvcfs
    #return {'gvcfs': gvcfs, 'gvcfs_idx': gvcfs_idx, 'doneFiles': doneFiles}

def make_intervals(output, outputDir, intDir, wildcards, dict_file, max_intervals):
    """Creates interval list files for parallelizing haplotypeCaller and friends. Writes one contig/chromosome per list file."""
    with open(dict_file, "r") as f:  # Read dict file to get contig info
        contigs = defaultdict()
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                contigs[chrom] = ln
        
        interval_file = output[0]
        
        with open(interval_file, "w") as fh:
            for contig, ln in contigs.items():
                print(f"{contig}\t1\t{ln}", file=fh)
        small_contigs = []
        not_small_contigs = []
        for i, (contig, ln) in enumerate(contigs.items()):
            if ln <= 500000:
                
                small_contigs.append((contig, ln))
            else:
                
                not_small_contigs.append((contig, ln))
        
        interval_list_file = Path(output.l, f"list0.list")
        with open(interval_list_file, "w") as fh:
            for contig, ln in small_contigs:
                print(f"{contig}:1-{ln}", file=fh)

        for i, (contig, ln) in enumerate(not_small_contigs):
            interval_list_file = Path(output.l, f"list{i+1}.list")
            with open(interval_list_file, "w") as fh:
                print(f"{contig}:1-{ln}", file=fh)
    
        '''
        This is the old code that was used to create the interval list files, with some max number of intervals desired. This will not be used for CCGP workflow.        
        else:
            ln_sum = sum(contigs.values())
            bp_per_interval = ln_sum // int(max_intervals)
            int_file = 0
            running_bp_total = 0
            out = deque()

            for chrom, ln in contigs.items():
                out.append(f"{chrom}:1-{ln}")
                running_bp_total += ln
                if running_bp_total >= bp_per_interval:
                    interval_file = Path(output.l, f"list{int_file}.list")
                    with open(interval_file, "a+") as f:
                        for _ in range(len(out)):
                            line = out.popleft()
                            print(line, file=f)
                    int_file += 1
                    running_bp_total = 0
            if out:
                interval_file = Path(output.l, f"list{int_file}.list")
                with open(interval_file, "a+") as f:
                    for _ in range(len(out)):
                        line = out.popleft()
                        print(line, file=f)
                int_file += 1
        '''
