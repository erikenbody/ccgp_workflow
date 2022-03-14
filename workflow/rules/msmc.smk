localrules: create_msmc_intervals

checkpoint create_msmc_intervals:
    input:
        bamcomplete = config['output'] + "{Organism}/{refGenome}/" + "bam_sumstats.tsv",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output: 
        directory(config['output'] + "{Organism}/{refGenome}/" + config["intDir"]),
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        LISTS = make_intervals_msmc(config["output"], config["intDir"], wildcards, input.dictf)

rule mappability:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed"
        #intdir = config['output'] + "{Organism}/{refGenome}/" + config["intDir"]
    conda:
        "../envs/msmc.yml"
    output: 
        fasta = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "L{list}.fasta",
        mapbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.bed",
        idx1 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.{ext}.{ext2}", ext=['ids', 'info', 'txt'], ext2=['concat', 'limits']),
        idx2 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        idx3 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.lf.{ext}.sbl", ext=['drp', 'drv']),
        idx4 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.rev.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
        idx5 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.rev.lf.{ext}.sbl", ext=['drp', 'drv']),
        idx6 = expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir'] + "index_{{list}}/" + "index.sa.{ext}", ext=['ind', 'len', 'val'])
    log:
        "logs/{Organism}/mappability/{refGenome}_{list}.txt"
    params:
        index = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['msmcDir'], "index_{list}"),
        bedprefix = os.path.join(workflow.default_remote_prefix, config['output'],  "{Organism}/{refGenome}/", config['msmcDir'], "L{list}"),
        l = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['intDir'], "list{list}.list"),
        ldir = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['intDir'])
    shell:
        """
        mkdir -p {params.ldir}
        gsutil cp {params.l} {params.ldir}
        samtools faidx {input.ref} -r {params.l} > {output.fasta}
        rm -r {params.index} && genmap index -F {output.fasta} -I {params.index}
        genmap map -K 150 -E 1 -I {params.index} -O {params.bedprefix} -b -T 1 #using -k 150 to approximate read length
        """

rule processmapp:
    input:
        fasta = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "L{list}.fasta",
        mapbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.bed",
    output:
        badbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.map.bed.gz"
    conda:
        "../envs/msmc.yml"
    log:
        "logs/{Organism}/mappability/{refGenome}_{list}.txt"
    shell:
        """
        awk '{{ if ($5 > 0.33 ) {{ print }} }}' {input.mapbed} | gzip > {output.badbed} 
        """

# rule genmap_index:
#     input:
#         ref = config["refGenomeDir"] + "{refGenome}.fna",
#     log:
#         "logs/{refGenome}/genmap_index.log"
#     conda:
#         "../envs/genmap.yml"
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
#     params:
#         os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap_index"))
#     output:
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.{ext}.{ext2}", ext=['ids', 'info', 'txt'], ext2=['concat', 'limits']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext}.sbl", ext=['drp', 'drv']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext}.sbl", ext=['drp', 'drv']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.sa.{ext}", ext=['ind', 'len', 'val'])
#     shell:
        # snakemake creates the output directory before the shell command, but genmap doesnt like this. so we remove the directory first.
#         "rm -rf {params} && genmap index -F {input.ref} -I {params} &> {log}"

# rule genmap_map:
#     input:
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.{ext}.{ext2}", ext=['ids', 'info', 'txt'], ext2=['concat', 'limits']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.lf.{ext}.sbl", ext=['drp', 'drv']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext2}", ext2=['drp', 'drs', 'drv', 'pst']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.rev.lf.{ext}.sbl", ext=['drp', 'drv']),
#         expand(config['output'] + "{{refGenome}}/" + "genmap_index/" + "index.sa.{ext}", ext=['ind', 'len', 'val'])
#     log:
#         "logs/{refGenome}/genmap_map.log"
#     params:
#         indir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap_index")),
#         outdir = os.path.join(workflow.default_remote_prefix, (config['output'] + "{refGenome}/" + "genmap"))
#     conda:
#         "../envs/genmap.yml"
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * res_config['genmap']['mem']
#     threads:
#         res_config['genmap']['threads']
#     output:
#         temp(config['output'] + "{refGenome}/" + "genmap/{refGenome}.genmap.bed")
#     shell:
#         """
#         genmap map -K 150 -E 0 -I {params.index} -O {params.bedprefix} -b -T {threads}  2>> {log} #using -k 150 to approximate read length
#         awk '{{ if ($5 > 0.33 ) {{ print }} }}' {params.bedprefix}.bed 2>> {log} | gzip > {output.badbed} 

#         """

rule msmc_input:
    """
    bcftools call for variant calling 1 individual 
    """
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        bam = lambda w: config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + samples['BioSample'].tolist()[0] + config['bam_suffix'], #this failed when I temped the final bam
        badbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.map.bed.gz",
        bamcaller = "msmctools/bamCaller.py",
        hetsep = "msmctools/generate_multihetsep.py"
    output:
        vcf = config['output'] + "{Organism}/{refGenome}/" + config["gvcfDir"]  + "{list}.vcf.gz",
        bed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "{list}.call.bed.gz",
        inmsc = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "{list}.msmc.input"
    params:
        l = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['intDir'], "list{list}.list")
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        SCAFF=$(cat {params.l})
        #DEPTH=$(samtools depth -r $SCAFF {input.bam} | awk '{{sum += $3}} END {{print sum / NR}}')
        DEPTH=2
        bcftools mpileup -q 20 -Q 20 -C 50 -r $SCAFF -f {input.ref} {input.bam} | bcftools call --ploidy 2 -c -V indels | {input.bamcaller} $DEPTH {output.bed} | gzip -c > {output.vcf}
        {input.hetsep} --mask={output.bed} --mask={input.badbed} {output.vcf} > {output.inmsc}
        """

rule gather_msmc_input:
    """
    #run the development version of msmc2
    """
    input:
        get_gather_msmc
    output: 
        fof = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc.fof",
        fof_clean = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc_clean.fof"
    conda:
        "../envs/msmc.yml"
    shell:
        """
        echo {input} > {output.fof}

        for i in `cat {output.fof}`
        do
            if [ -s $i ]; then
                realpath $i >> {output.fof_clean}
            fi
        done
        """

# rule run_msmc:
#     input: 
#         get_gather_msmc
#     output: 
#         mscmc = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc.final.txt"
#     params:
#         prefix = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc"
#     conda:
#         "../envs/msmc.yml"
#     threads: 
#         res_config['run_msmc2']['threads']
#     benchmark:
#         "benchmarks/{Organism}/sortVcf/{refGenome}.txt"
#     shell:
#         """
#         MTOOLS=/scratch/home/eenbody/tools/msmc-tools
#         $MTOOLS/msmc_1.1.0_linux64bit -t {threads} --outFilePrefix {params.prefix} {input}
#         """

rule run_msmc2:
    """
    #run the development version of msmc2
    """
    input: 
        fof_clean = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc_clean.fof",
        msmctools = "msmctools/msmc2_linux64bit"
    output: 
        mscmc = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.txt",
        mscmcD = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.default.final.txt"
    params:
        prefix2 = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", "{Organism}_{refGenome}_msmc2"),
        prefixD = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", "{Organism}_{refGenome}_msmc2.default")
    conda:
        "../envs/msmc.yml"
    threads:
        res_config['run_msmc2']['threads']
    shell:
        """
        chmod +x {input.msmctools}
        {input.msmctools} -t {threads} -p 1*4+25*2+1*4+1*6 --outFilePrefix {params.prefix2} $(<{input.fof_clean})
        {input.msmctools} -t {threads} --outFilePrefix {params.prefixD} $(<{input.fof_clean})
        """

rule msmc_plots:
    """
    Call plotting script
    """
    input:
        mscmc = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.txt"
    output: 
        mscmcPlot = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.pdf"
    conda:
        "../envs/msmc.yml"
    script:
        "../scripts/plot_msmc.R"


#cleanup script is problematic, because if bam is input, the rerunning the pipeline (e.g. editing plotting script) leads to wanting to re run the full pipeline
# rule cleanup:
#     """
#     delete bam file
#     """
#     input:
#         bam = lambda w: config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + samples['BioSample'].tolist()[0] + config['bam_suffix'], #this failed when I temped the final bam
#         mscmcPlot = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.pdf"
#     output:
#         clean = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_cleanup.txt"
#     shell:
#         """
#         if [ -s {input.mscmcPlot} ]; then
#             rm {input.bam}
#         fi
#         echo "cleanup success" > {output.clean}
#         """