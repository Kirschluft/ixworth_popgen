rule trim_reads:
    input:
        unpack(get_fastq),
    output:
        trimmed=[
            temp("results/trimmed/{sample}-{unit}.1.fastq.gz"),
            temp("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        ],
        unpaired1=temp("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        unpaired2=temp("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        html="results/qc/fastp/{sample}-{unit}_pe.html",
        json="results/qc/fastp/{sample}-{unit}_pe.json",
    params:
        extra=config["params"]["fastp"] if config["params"]["fastp"] else "",
    threads: 8
    log:
        "logs/fastp/{sample}-{unit}.log",
    conda:
        "../envs/fastp.yaml"
    wrapper:
        "v5.1.0/bio/fastp"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}-{unit}.bam"),
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="none",
    threads: 8
    wrapper:
        "v5.1.0/bio/bwa/mem"


rule merge_bam:
    input:
        lambda wildcards: [
            f"results/mapped/{wildcards.sample}-{unit}.bam"
            for unit in get_units(wildcards)
        ],
    output:
        temp("results/merged/{sample}.bam"),
    log:
        "logs/samtools_merge/{sample}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v5.1.0/bio/samtools/merge"


rule sort_bam:
    input:
        "results/merged/{sample}.bam",
    output:
        temp("results/merged/{sample}.sorted.bam"),
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        """
          samtools sort -o {output} {input} > {log} 2>&1
          samtools index {output} > {log} 2>&1
        """


rule mark_duplicates:
    input:
        bams="results/merged/{sample}.sorted.bam",
    output:
        bam=temp("results/dedup/{sample}.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        f'REMOVE_DUPLICATES={config["processing"]["remove-duplicates"]}',
    resources:
        mem_mb=16000,
    wrapper:
        "v5.1.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/chromosomes.fasta",
        dict="resources/chromosomes.dict",
        known="resources/variation.vcf.gz",
        known_idx="resources/variation.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp"),
    log:
        "logs/gatk/bqsr/{sample}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    resources:
        mem_mb=16000,
    wrapper:
        "v5.1.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/chromosomes.fasta",
        dict="resources/chromosomes.dict",
        recal_table="results/recal/{sample}.grp",
    output:
        bam=temp("results/recal/{sample}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}.log",
    params:
        extra=get_regions_param(),
    resources:
        mem_mb=16000,
    wrapper:
        "v5.1.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        temp("{prefix}.bam.bai"),
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "v5.1.0/bio/samtools/index"
