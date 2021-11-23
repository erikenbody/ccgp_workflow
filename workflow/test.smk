rule all:
    input: "results/limpet/xgLotScab1.0.p_ctg/01_mappedReads/M0D059136X_final.bam"
    output: "test/b.txt"
    benchmark:
        "test/benchmark/bb.txt"
    log:
        "test/log/b.txt"
    shell:
        "mv {input} {output}"