# workaround since fastqc can not handle paired end reads
rule fastqc_read1:
    input:
        unpack(get_read1),
    output:
        html="results/qc/fastqc/{sample}-{unit}-1_fastqc.html",
        zip="results/qc/fastqc/{sample}-{unit}-1_fastqc.zip",
    log:
        "logs/fastqc/{sample}-{unit}-1.log",
    resources:
        mem_mb=10000,
    threads: 2
    wrapper:
        "v5.1.0/bio/fastqc"


rule fastqc_read2:
    input:
        unpack(get_read2),
    output:
        html="results/qc/fastqc/{sample}-{unit}-2_fastqc.html",
        zip="results/qc/fastqc/{sample}-{unit}-2_fastqc.zip",
    log:
        "logs/fastqc/{sample}-{unit}-2.log",
    resources:
        mem_mb=10000,
    threads: 2
    wrapper:
        "v5.1.0/bio/fastqc"


rule qualimap_stats_wgs:
    input:
        bam="results/dedup/{sample}.bam",
    output:
        directory("results/qc/qualimap/{sample}"),
    log:
        "logs/bamqc/{sample}.log",
    threads: 4
    resources:
        mem_mb=16000,
    wrapper:
        "v5.1.0/bio/qualimap/bamqc"


rule samtools_stats:
    input:
        "results/dedup/{sample}.bam",
    output:
        "results/qc/samtools-stats/{sample}.txt",
    log:
        "logs/samtools-stats/{sample}.log",
    wrapper:
        "v5.1.0/bio/samtools/stats"


rule samtools_idxstats:
    input:
        bam="results/dedup/{sample}.bam",
        bai="results/dedup/{sample}.bam.bai",
    output:
        "results/qc/samtools-idxstats/{sample}.idxstats.txt",
    log:
        "logs/samtools-idxstats/{sample}.log",
    wrapper:
        "v5.1.0/bio/samtools/idxstats"


rule samtools_depth:
    input:
        bams=lambda wildcards: [f"results/dedup/{wildcards.sample}.bam"],
    output:
        "results/qc/samtools-depth/{sample}.depth.txt",
    log:
        "logs/samtools-depth/{sample}.log",
    wrapper:
        "v5.1.0/bio/samtools/depth"


rule multiqc:
    input:
        expand(
            "results/qc/fastqc/{u.sample}-{u.unit}-1_fastqc.zip", u=units.itertuples()
        ),
        expand(
            "results/qc/fastqc/{u.sample}-{u.unit}-2_fastqc.zip", u=units.itertuples()
        ),
        expand("results/qc/fastp/{u.sample}-{u.unit}_pe.html", u=units.itertuples()),
        expand("results/qc/samtools-stats/{u.sample}.txt", u=samples.itertuples()),
        expand(
            "results/qc/samtools-idxstats/{u.sample}.idxstats.txt",
            u=samples.itertuples(),
        ),
        expand(
            "results/qc/samtools-depth/{u.sample}.depth.txt",
            u=samples.itertuples(),
        ),
        expand("results/qc/qualimap/{u.sample}", u=samples.itertuples()),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    wrapper:
        "v5.1.0/bio/multiqc"
