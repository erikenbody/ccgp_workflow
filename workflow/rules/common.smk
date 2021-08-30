import glob
import re
import os

### INPUT FUNCTIONS ###
def get_fastp_stats(wildcards):
    checkpoint_output = checkpoints.get_fastq_pe.get(**wildcards).output[0]
    fastq_dir = os.path.join(config["fastqDir"], wildcards.Organism, wildcards.run, "*.fastq")
    fastq_files = glob.glob(fastq_dir)

    if len(fastq_files) == 2:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{run}/{run}.out", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())
    else:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['sumstatDir'] + "{run}/se_{run}.out", run=samples.loc[samples['BioSample'] == wildcards.sample]['Run'].tolist())

def get_reads(wildcards):
    checkpoint_output = checkpoints.get_fastq_pe.get(**wildcards).output[0]
    fastq_dir = os.path.join(config["fastqDir"], wildcards.Organism, wildcards.run, "*.fastq")
    fastq_files = glob.glob(fastq_dir)

    print(fastq_dir, fastq_files)
    if len(fastq_files) == 1:
        return config['output'] + "{Organism}/{refGenome}/" + config['fastqFilterDir'] + "{sample}/{run}/{run}.fastq.gz"
    else:
        return expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['fastqFilterDir'] + "{{sample}}/{{run}}/{{run}}_{num}.fastq.gz", num=[1,2])

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
    gvcfs = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz", sample=sample_names)
    gvcfs_idx = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.raw.g.vcf.gz.tbi", sample=sample_names)
    doneFiles = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['gvcfDir'] + "{sample}/" + "L{{list}}.done", sample=sample_names)
    return {'gvcfs': gvcfs, 'gvcfs_idx': gvcfs_idx, 'doneFiles': doneFiles}

def make_intervals(outputDir, intDir, wildcards, dict_file, num_splits):
    listfile_index = 0
    with open(dict_file, "r") as f:
        for line in f:
            if line.startswith("@SQ"):
                line = line.split("\t")
                chrom = line[1].split("SN:")[1]
                ln = int(line[2].split("LN:")[1])
                step_size = ln // int(num_splits)
                nums = list(range(1, ln, step_size))
                lists = [[nums[i-1], nums[i]] for i in range(1, len(nums))]
                lists[-1][1] = ln
                for l in lists:
                    if l[0] != 1:
                        l[0] += 1
                    list_file = gatk_list_dir = os.path.join(outputDir, wildcards.Organism, wildcards.refGenome, intDir, f"list{listfile_index}.list")
                    with open(list_file, "w") as f:
                        print(f"{chrom}:{l[0]}-{l[1]}", file=f)
                    listfile_index += 1
                