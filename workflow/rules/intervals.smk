localrules: picard_intervals, create_intervals

rule picard_intervals:
    input:
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
    output:
        intervals = temp(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list")
    params:
        minNmer = int(config['minNmer'])
    conda:
        '../envs/intervals.yml'
    log:
        "logs/{Organism}/{refGenome}/picard_intervals/log"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['process_ref']['mem']   
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} &> {log}\n" 

checkpoint create_intervals:
    input:
        fai = config["refGenomeDir"] + "{refGenome}.fna" + ".fai",
        dictf = config["refGenomeDir"] + "{refGenome}" + ".dict",
        ref = config["refGenomeDir"] + "{refGenome}.fna",
        intervals = ancient(config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_output.interval_list")
    params:
        maxIntervalLen = int(config['maxIntervalLen']),
        maxBpPerList = int(config['maxBpPerList']),
        maxIntervalsPerList = int(config['maxIntervalsPerList']),
        minNmer = int(config['minNmer']),
        max_intervals = config['maxNumIntervals']
    output: 
        config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "{refGenome}_intervals_fb.bed",
        #config['output'] + "{Organism}/{refGenome}/" + config["intDir"] + "num_intervals.txt",
        l = directory(config['output'] + "{Organism}/{refGenome}/" + config['intDir'] + "lists/")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * res_config['create_intervals']['mem'] 
    run:
        Path(output.l).mkdir(exist_ok=True)
        if config['split_by_n']:
            NUM_OF_LISTS = helperFun.createListsGetIndices(params.maxIntervalLen, params.maxBpPerList, params.maxIntervalsPerList, params.minNmer, config["output"], config["intDir"], wildcards, input.dictf, input.intervals)
        else:
            NUM_OF_LISTS = make_intervals(output, config["output"], config["intDir"], wildcards, input.dictf, params.max_intervals)
