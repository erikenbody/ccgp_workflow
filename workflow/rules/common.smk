import glob
import re
import os
import pandas as pd
from collections import defaultdict, deque
from snakemake.exceptions import WorkflowError
### INPUT FUNCTIONS ###

def write_coords_file(wildcards):
    out_df = samples[["BioSample", "long", "lat"]]
    out_df.drop_duplicates("BioSample", inplace=True)
    outpath = os.path.join(workflow.default_remote_prefix, config['output'], f"{wildcards.Organism}/{wildcards.refGenome}/", config['qcDir'], f"{wildcards.Organism}_{wildcards.refGenome}.coords.txt")
    out_df.to_csv(outpath, index=False, sep="\t")
def get_coords_if_available(wildcards):
    if 'lat' in samples.columns and 'long' in samples.columns:
        return config['output'] + "{Organism}/{refGenome}/" + config['qcDir'] + "{Organism}_{refGenome}.coords.txt" 
def get_ena_url(wildcards):
    prefix = wildcards.run[:6]
    lastdigit = wildcards.run[-1]
    code = wildcards.run[:3]
    baseloc = "http://ftp.sra.ebi.ac.uk/vol1/"
    if len(wildcards.run) > 9:
        url = prefix + "/" + "00" + lastdigit + "/" + wildcards.run
    else:
        url = prefix + "/" + wildcards.run
    sra_url = baseloc + code + "/" + url
    fastq_url = baseloc + "fastq/" + url
    return {"sra_url": sra_url, "fastq_url": fastq_url}

def get_gvcf_cmd(wildcards):
    sample_names = samples.BioSample.unique().tolist()
    vcfs = expand(os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['gvcfDir'], "{sample}.g.vcf.gz"), **wildcards, sample=sample_names)
    out = " ".join(["-v " + vcf for vcf in vcfs])
    return out
def get_tbis_for_list(wildcards):
    sample_names = samples.BioSample.unique().tolist()
    return expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz.tbi", **wildcards, sample=sample_names)

def get_gvcfs_for_list(wildcards):
    sample_names = samples.BioSample.unique().tolist()
    return expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}/" + "L{list}.raw.g.vcf.gz", **wildcards, sample=sample_names)

def get_gvcfs(wildcards):
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].unique().tolist()
    return expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz", **wildcards, sample=_samples)

def get_tbis(wildcards):
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].unique().tolist()
    return expand(config['output'] + "{Organism}/{refGenome}/" + config['gvcfDir'] + "{sample}.g.vcf.gz.tbi", **wildcards, sample=_samples)

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

def get_bams_for_dedup(wildcards):
    runs = samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist()
    if len(runs) == 1:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{{sample}}/{run}.bam", run=runs)
    else:
        return config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{sample}.bam"

def get_reads(wildcards):
    """Returns local read files if present. Defaults to SRR if no local reads in sample sheet."""
    if config['remote_reads']:
        return get_remote_reads(wildcards)
    else:
        row = samples.loc[samples['Run'] == wildcards.run]
        if 'fq1' in samples.columns and 'fq2' in samples.columns:
            if os.path.exists(row.fq1.item()) and os.path.exists(row.fq2.item()):
                r1 = row.fq1.item()
                r2 = row.fq2.item()
                return {"r1": r1, "r2": r2}
            else:
                raise WorkflowError(f"fq1 and fq2 specified for {wildcards.sample}, but files were not found.")
        else:
            r1 = config["fastqDir"] + f"{wildcards.Organism}/{wildcards.sample}/{wildcards.run}_1.fastq.gz",
            r2 = config["fastqDir"] + f"{wildcards.Organism}/{wildcards.sample}/{wildcards.run}_2.fastq.gz"
        return {"r1": r1, "r2": r2}

def get_remote_reads(wildcards):
    """Use this for reads on a different remote bucket than the default."""
    row = samples.loc[samples['Run'] == wildcards.run]
    r1 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq1.item()))
    r2 = GS.remote(os.path.join(GS_READS_PREFIX, row.fq2.item()))
    return {"r1": r1, "r2": r2}



def get_read_group_sentieon(wildcards):
    """Denote sample name and library_id in read group."""
    return r"'@RG\tID:{lib}\tSM:{sample}\tPL:ILLUMINA'".format(
        sample=wildcards.sample,
        lib=wildcards.sample
    )

def get_read_group_bwa(wildcards):
    """Denote sample name and library_id in read group."""
    return r"-R '@RG\tID:{lib}\tSM:{sample}\tPL:ILLUMINA'".format(
        sample=wildcards.sample,
        lib=samples.loc[samples['BioSample'] == wildcards.sample]["LibraryName"].tolist()[0]
    )

def check_contig_names(fai, touch_file):

    dffai = pd.read_table(fai, sep='\t', header = None)
    fai_result=pd.to_numeric(dffai[0], errors='coerce').notnull().all()
    if fai_result==True:
        print("QC plots not generated because contig names are numeric and plink does not accept numeric contig names")
    elif fai_result==False:
        with open(touch_file, "w") as writer:
            writer.write("contigs are strings")

def get_sumstats(wildcards):
    # Gets the correct sample given the organism and reference genome
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    fastpFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_fastp.out", sample=_samples)
    alnSumMetsFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_AlnSumMets.txt", sample=_samples)
    coverageFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_coverage.txt", sample=_samples)
    return {'alnSumMetsFiles': alnSumMetsFiles, 'coverageFiles': coverageFiles, 'fastpFiles': fastpFiles}

def get_aln_sum_metrics(wildcards):
    """Snakemake seems to struggle with unpack() + expand() and default_remote_prefix. So we have to do this one by one.
       See snakemake/snakemake #1396 for details."""
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_AlnSumMets.txt", sample=_samples)

def get_fastpFiles(wildcards):
    """Snakemake seems to struggle with unpack() + expand() and default_remote_prefix. So we have to do this one by one.
       See snakemake/snakemake #1396 for details."""
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_fastp.out", sample=_samples)

def get_coverageFiles(wildcards):
    """Snakemake seems to struggle with unpack() + expand() and default_remote_prefix. So we have to do this one by one.
       See snakemake/snakemake #1396 for details."""
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{sample}_coverage.txt", sample=_samples)

def get_gather_vcfs(wildcards):
    """
    Gets filtered vcfs for gathering step. This function gets the interval list indicies from the corresponding
    genome, then produces the file names for the filtered vcf with list index."""
    checkpoint_output = checkpoints.create_intervals.get(**wildcards).output[0]
    list_dir_search = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['intDir'], "*.list")
    list_files = glob.glob(list_dir_search)
    out = []
    for f in list_files:
        f = os.path.basename(f)
        index = re.search("\d+", f).group() # Grab digits from list file name and put in out list
        vcf = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['vcfDir_gatk'], f"filtered_L{index}.vcf")
        out.append(vcf)

    return out

def gather_vcfs_CLI(wildcards):
    """
    Gatk enforces that you have a -I before each input vcf, so this function makes that string
    """
    vcfs = get_gather_vcfs(wildcards)
    out = " ".join(["-I " + vcf for vcf in vcfs])
    return out

def write_db_mapfile(wildcards):
    dbMapFile = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['dbDir'], f"DB_mapfile_L{wildcards.list}")
    sample_names = set(samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist())
    with open(dbMapFile, 'w') as f:
        for sample in sample_names:
            gvcf_path = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['gvcfDir'], sample, f"L{wildcards.list}.raw.g.vcf.gz")
            print(sample, gvcf_path, sep="\t", file=f)

def get_input_for_mapfile(wildcards):
    sample_names = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].tolist()
    gvcfs = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz", sample=sample_names, **wildcards)
    gvcfs_idx = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz.tbi", sample=sample_names, **wildcards)
    doneFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.done", sample=sample_names, **wildcards)
    # only gvcfs are needed for the mapfile
    return gvcfs

def get_bedgraph_to_convert(wildcards):
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].unique().tolist()
    if len(_samples) == 1:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{sample}.sorted.bg", sample=_samples)
    else:
        return config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + "postMerge/{Organism}.merge.bg"

def get_input_for_coverage(wildcards):
    # Gets the correct sample given the organism and reference genome for the bedgraph merge step
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].unique().tolist()
    bedgraphFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{sample}" + ".sorted.bg", sample=_samples)
    chromFile = config['output'] + "{refGenome}/" + "{refGenome}" + ".sizes"
    return {'bedgraphs': bedgraphFiles, 'chrom': chromFile}

def get_bedgraphs(wildcards):
    """Snakemake seems to struggle with unpack() and default_remote_prefix. So we have to do this one by one."""
    _samples = samples.loc[(samples['Organism'] == wildcards.Organism) & (samples['refGenome'] == wildcards.refGenome)]['BioSample'].unique().tolist()
    bedgraphFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['bamDir'] + "preMerge/{sample}" + ".sorted.bg", sample=_samples)
    return bedgraphFiles

def make_intervals(outputDir, intDir, wildcards, dict_file, max_intervals):
    """Creates interval list files for parallelizing haplotypeCaller and friends. Writes one contig/chromosome per list file."""

    with open(dict_file, "r") as f:  # Read dict file to get contig info
        contigs = defaultdict()
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                contigs[chrom] = ln

        interval_file = os.path.join(outputDir,wildcards.Organism,wildcards.refGenome,intDir, f"{wildcards.refGenome}_intervals_fb.bed")
        with open(interval_file, "w") as fh:
            for contig, ln in contigs.items():
                print(f"{contig}\t1\t{ln}", file=fh)

        if len(contigs.values()) <= max_intervals:
            for i, (contig, ln) in enumerate(contigs.items()):
                interval_list_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{i}.list")
                with open(interval_list_file, "w") as f:
                    print(f"{contig}:1-{ln}", file=f)

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
                    interval_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{int_file}.list")
                    with open(interval_file, "a+") as f:
                        for _ in range(len(out)):
                            line = out.popleft()
                            print(line, file=f)
                    int_file += 1
                    running_bp_total = 0
            if out:
                interval_file = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{int_file}.list")
                with open(interval_file, "a+") as f:
                    for _ in range(len(out)):
                        line = out.popleft()
                        print(line, file=f)

def make_intervals_msmc(outputDir, intDir, wildcards, dict_file):
    """Creates interval list files for parallelizing haplotypeCaller and friends. Writes one contig/chromosome per list file."""
    
    with open(dict_file, "r") as f:  # Read dict file to get contig info
        contigs = defaultdict()
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                contigs[chrom] = ln
        
        interval_file = os.path.join(workflow.default_remote_prefix,outputDir,wildcards.Organism,wildcards.refGenome,intDir, f"{wildcards.refGenome}_msmc_intervals_fb.bed")
        with open(interval_file, "w") as fh:
            for contig, ln in contigs.items():
                print(f"{contig}\t1\t{ln}", file=fh)
        #this is a little ugly because "i" includes the full unfiltered number of scaffolds
        for i, (contig, ln) in enumerate(contigs.items()):
            if ln > 1:
                interval_list_file = os.path.join(workflow.default_remote_prefix, outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{i}.list")
                with open(interval_list_file, "w") as f:
                    print(f"{contig}:1-{ln}", file=f)

def get_gather_msmc(wildcards):
    """
    Gets msmc formatted files for gathering step. This function gets the interval list indicies from the corresponding
    genome, then produces the file names for the filtered vcf with list index."""
    checkpoint_output = checkpoints.create_msmc_intervals.get(**wildcards).output[0]
    list_dir_search = os.path.join(workflow.default_remote_prefix, config['output'], wildcards.Organism, wildcards.refGenome, config['intDir'], "*.list")
    list_files = glob.glob(list_dir_search)
    out = []
    for f in list_files:
        f = os.path.basename(f)
        index = re.search("\d+", f).group() # Grab digits from list file name and put in out list
        msmc = os.path.join(config['output'], wildcards.Organism, wildcards.refGenome, config['msmcDir'], f"{index}.msmc.input")
        out.append(msmc)
    return out
