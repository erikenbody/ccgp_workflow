localrules: create_msmc_intervals, gather_msmc_input

rule create_msmc_intervals:
    input:
        bamcomplete = config['output'] + "{Organism}/{refGenome}/" + "bam_sumstats.tsv",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna"
    output:
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed",
        expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config["intDir"] + "list{list}.list", list=range(config['maxNumIntervals']))
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        LISTS = make_intervals(config["output"], config["intDir"], wildcards, input.dictf, config['maxNumIntervals'])

rule mappability:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed",
        list_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "list{list}.list"
    conda:
        "../envs/msmc.yml"
    output: 
        fasta = config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "L{list}.fasta",
        mapbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.bed"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['mappability']['mem']
    log:
        "logs/{Organism}/mappability/{refGenome}_{list}.txt"
    params:
        index = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['msmcDir'], "index_{list}"),
        bedprefix = os.path.join(workflow.default_remote_prefix, config['output'],  "{Organism}/{refGenome}/", config['msmcDir'], "L{list}"),
        l = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['intDir'], "list{list}.list"),
        ldir = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['intDir'])
    shell:
        """
        bedtools getfasta -fi {input.ref} -bed {input.list_file} -fo {output.fasta}
        rm -rf {params.index} && genmap index -F {output.fasta} -I {params.index}
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

rule bamdepth:
    input:
        bam = lambda w: config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + samples['BioSample'].tolist()[0] + config['bam_suffix'], #this failed when I temped the final bam
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed",
        list_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "list{list}.list"
    output:
        depth = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "bam{list}_depth.txt"
    conda:
        "../envs/msmc.yml"
    shell:
        """
        set +u
        samtools depth -b {input.list_file} {input.bam} > {output.depth}
        #DEPTH=$(awk "{{sum += \$3}} END {{print sum / NR}}" {output.depth})
        #>&2 echo $DEPTH
        """

rule msmc_input:
    """
    bcftools call for variant calling 1 individual 
    """
    input:
        depth = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "bam{list}_depth.txt",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        bam = lambda w: config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + samples['BioSample'].tolist()[0] + config['bam_suffix'], #this failed when I temped the final bam
        bai = lambda w: config['output'] + "{Organism}/{refGenome}/" + config['bamDir'] + samples['BioSample'].tolist()[0] + config['bam_suffix'] + '.bai', #this failed when I temped the final bam
        badbed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir'] + "L{list}.map.bed.gz",
        msmctools = get_msmc_tools,
        intervals = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_msmc_intervals_fb.bed",
        list_file = config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "list{list}.list"
    output:
        vcf = config['output'] + "{Organism}/{refGenome}/" + config["gvcfDir"]  + "{list}.vcf.gz",
        bed = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "{list}.call.bed.gz",
        inmsc = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "{list}.msmc.input",
        mdir = directory(config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "{list}.mscmc_input_files")
    params:
        bamcaller = os.path.join(workflow.default_remote_prefix, "msmctools/bamCaller.py",),
        hetsep = os.path.join(workflow.default_remote_prefix, "msmctools/generate_multihetsep.py",),
        prefix = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['msmcDir'], "{list}.mscmc_input_files")
    conda:
        "../envs/bcftools.yml"
    resources:
        machine_type="n2-standard",
        mem_mb = lambda wildcards, attempt: attempt * res_config['msmc_input']['mem']
    shell:
        """
        set +u
        chmod +x {params.bamcaller}
        chmod +x {params.hetsep}
        chmod +x workflow/scripts/mean_depth.py

        DEPTH=$(workflow/scripts/mean_depth.py {input.depth})
        >&2 echo $DEPTH

        mkdir -p {params.prefix}

        while read chr start end
        do
            bcftools mpileup -q 20 -Q 20 -C 50 -r ${{chr}}:${{start}}-${{end}} -f {input.ref} {input.bam} | bcftools call --ploidy 2 -c -V indels | {params.bamcaller} $DEPTH {output.bed} | gzip -c > {output.vcf}
            {params.hetsep} --mask={output.bed} --mask={input.badbed} {output.vcf} > {params.prefix}/msmc_input_${{chr}}.input

            if [ -s {params.prefix}/msmc_input_${{chr}}.input ]; then
                echo {params.prefix}/msmc_input_${{chr}}.input >> {output.inmsc}
            fi

        done < {input.list_file}

        """

rule gather_msmc_input:
    input:
        get_gather_msmc
    output: 
        #fof = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc.fof",
        fof_clean = config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc_clean.fof"
    conda:
        "../envs/msmc.yml"
    shell:
        """
        cat {input} > {output.fof_clean}
        """

rule run_msmc2:
    """
    #run the development version of msmc2
    """
    input:
        fof_clean = ancient(config['output'] + "{Organism}/{refGenome}/" + config['msmcDir']  + "msmc_clean.fof"),
        msmctools = get_msmc_tools,
        #clean_files = get_clean_files,
        mdir = directory(expand(config['output'] + "{{Organism}}/{{refGenome}}/" + config['msmcDir']  + "{list}.mscmc_input_files", list=range(config["maxNumIntervals"])))
    output: 
        mscmc = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.final.txt",
        mscmcD = config['output'] + "{Organism}/{refGenome}/" + "{Organism}_{refGenome}_msmc2.default.final.txt"
    params:
        # direct = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['msmcDir']),
        prefix2 = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", "{Organism}_{refGenome}_msmc2"),
        prefixD = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", "{Organism}_{refGenome}_msmc2.default"),
        tool = os.path.join(workflow.default_remote_prefix,"msmctools/msmc2_linux64bit"),
        #inpath = os.path.join(workflow.default_remote_prefix, config['output'], "{Organism}/{refGenome}/", config['msmcDir'], "mscmc_input_files")
    conda:
        "../envs/msmc.yml"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['run_msmc2']['mem']
    threads:
        res_config['run_msmc2']['threads']
    shell:
        """
        chmod +x {params.tool}
        {params.tool} -t {threads} -p 1*4+25*2+1*4+1*6 --outFilePrefix {params.prefix2} $(<{input.fof_clean})
        {params.tool} -t {threads} --outFilePrefix {params.prefixD} $(<{input.fof_clean})
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
