rule intervals:
    """
    Because snakemake input MUST be a file, this hack is necessary to iterate through intervals in the next step,
    inputing one interval in the command at a time. This information is stored in intervals_fb
    """
    output:
        interval = temp(touch(intervalDir + "{i}"))

rule bam2vcf:
    """
    Call freebayes.
    """
    input:
        ref = config['ref'],
        bams = expand(bamDir + "{sample}" + bam_suffix, sample=SAMPLES) # input all bams at once
    output: 
        vcf = temp(vcfDir_fb + "{i}.vcf") # output one vcf per interval
    conda:
        "../envs/bam2vcf.yml"
    resources:
        # increase memory every attempt
        mem_mb = lambda wildcards, attempt: attempt * res_config['bam2vcf']['mem'] 
    shell:
        "freebayes -f {input.ref} "
        "-r {wildcards.i} "
        "-g {maxDP_fb} "
        "{input.bams} "
        "> {output.vcf}\n"

rule filterVcfs:
    """
    This rule filters all of the VCFs
    """
    input:
        vcf = vcfDir + "L{list}.vcf",
        ref = config['ref']
    output: 
        vcf = vcfDir + "L{list}_filter.vcf"
    params:
        gatherVcfsInput = helperFun.getVcfs_gatk(LISTS, vcfDir)
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['filterVcfs']['mem']   # this is the overall memory requested
    shell:
        "gatk VariantFiltration "
        "-R {input.ref} "
        "-V {input.vcf} " 
        "--output {output.vcf} "
        "--filter-name \"RPRS_filter\" "
        "--filter-expression \"(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)) || ((vc.isIndel() || vc.isMixed()) && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)) || (vc.hasAttribute('QD') && QD < 2.0)\" "
        "--filter-name \"FS_SOR_filter\" "
        "--filter-expression \"(vc.isSNP() && ((vc.hasAttribute('FS') && FS > 60.0) || (vc.hasAttribute('SOR') &&  SOR > 3.0))) || ((vc.isIndel() || vc.isMixed()) && ((vc.hasAttribute('FS') && FS > 200.0) || (vc.hasAttribute('SOR') &&  SOR > 10.0)))\" "
        "--filter-name \"MQ_filter\" "
        "--filter-expression \"vc.isSNP() && ((vc.hasAttribute('MQ') && MQ < 40.0) || (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))\" "
        "--filter-name \"QUAL_filter\" "
        "--filter-expression \"QUAL < 30.0\" "
        "--invalidate-previous-filters true\n"

rule gatherVcfs:
    """
    This rule gathers all the filtered VCFs, one per list, into one final VCF 
    """
    input:
        vcfs = expand(vcfDir + "L{list}_filter.vcf", list=LISTS)
    output:
        vcf = config["gatkDir"] + config["spp"] + ".vcf.gz"
    params:
        gatherVcfsInput = helperFun.getVcfs_gatk(LISTS, vcfDir)
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['gatherVcfs']['mem']   # this is the overall memory requested
    shell:
        "gatk GatherVcfs "
        "{params.gatherVcfsInput} "
        "-O {output.vcf}"

rule sortVcf:
    """
    Sort VCF for more efficient processing of vcftools and bedtools
    """
    input:
        vcf = config["gatkDir"] + config["spp"] + ".vcf.gz"
    output:
        vcf = config["gatkDir"] + config["spp"] + "_final.vcf.gz"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['sortVcf']['mem']   # this is the overall memory requested
    shell:
        "picard SortVcf -I {input.vcf} -O {output.vcf}"

rule vcftools:
    input:
        vcf = config["gatkDir"] + config['spp'] + "_final.vcf.gz",
        int = intDir + config["genome"] + "_intervals_fb.bed"
    output: 
        missing = gatkDir + "missing_data_per_ind.txt"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['vcftools']['mem']    # this is the overall memory requested
    shell:
        "vcftools --gzvcf {input.vcf} --remove-filtered-all --minDP 1 --stdout --missing-indv > {output.missing}"

rule bedtools:
    input:
        vcf = config["gatkDir"] + config['spp'] + "_final.vcf.gz",
        int = intDir + config["genome"] + "_intervals_fb.bed"
    output: 
        SNPsPerInt = gatkDir + "SNP_per_interval.txt"
    conda:
        "../envs/bam2vcf.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * res_config['bedtools']['mem']    # this is the overall memory requested
    shell:
        "bedtools coverage -a {input.int} -b {input.vcf} -counts -sorted > {output.SNPsPerInt}"
