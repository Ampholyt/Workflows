import os
NAMES = ["ERR024604_tiny","ERR024605_tiny",
        "ERR024606_tiny","ERR024607_tiny",
        "ERR024608_tiny","ERR024609_tiny"]

# dir = "/home/saloberg/Schreibtisch/Workflows/EinführungsWoche/SWP_Workflows/Tag1"
# NAMES, = glob_wildcards(os.path.join(dir,"/fastq/tiny/{dataset}_1.fastq.gz"))
# IDS, = glob_wildcards("fastq/tiny/{dataset}_{id}.fastq")

# REF = ["ref-genom.1","ref-genom.2","ref-genom.3","ref-genom.4","ref-genom.rev.1","ref-genom.rev.2"]
rule all:
    input:
        # "stats/{dataset}_stats.txt".format(dataset=NAMES[0])
        # "mapped/{dataset}.sam".format(dataset=NAMES[0])
        # ["mapped/{dataset}.sam".format(dataset=dataset) for dataset in NAMES]
        ["stats/tiny/{dataset}_stats.txt".format(dataset=dataset) for dataset in NAMES]
        # expand("stats/tiny/{dataset}_stats.txt", dataset=NAMES)

# rule ref_index:
#     input:
#         "reference/reference.fa"
#     output:
#         "ref-genom.1.bt2"
#         "ref-genom.2.bt2"
#         "ref-genom.3.bt2"
#         "ref-genom.4.bt2"
#         "ref-genom.rev.1.bt2"
#         "ref-genom.rev.2.bt2"
#     shell:
#         "bowtie2-build {input} reference/ref-genom

rule mapping:
    input:
        "fastq/tiny/{dataset}_1.fastq.gz",
        "fastq/tiny/{dataset}_2.fastq.gz"
    output:
        "mapped/tiny/{dataset}.sam"
    shell:
        "bowtie2 --threads 4 -x reference/ref-genom -1 {input[0]} -2 {input[1]} > {output}"


rule samToBam:
    input:
        "mapped/tiny/{dataset}.sam"
    output:
        "bam/tiny/{dataset}.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort:
    input:
        "bam/tiny/{dataset}.bam"
    output:
        "bam_sorted/tiny/{dataset}_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule indexing:
    input:
        "bam_sorted/tiny/{dataset}_sorted.bam"
    output:
        "bam_sorted/tiny/{dataset}_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule stats:
    input:
        "bam_sorted/tiny/{dataset}_sorted.bam",
        "bam_sorted/tiny/{dataset}_sorted.bam.bai"
    output:
        "stats/tiny/{dataset}_stats.txt"
    shell:
        "samtools idxstats {input[0]} > {output}"
