if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output} 2> {log}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/chromosomes.fasta",
        idx="resources/chromosomes.dict",
        known="resources/variation.vcf.gz",
        tbi="resources/variation.vcf.gz.tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=temp("results/called/{sample}.{contig}.g.vcf.gz"),
        tbi=temp("results/called/{sample}.{contig}.g.vcf.gz.tbi"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
    threads: 2
    resources:
        mem_mb=64000,
    wrapper:
        "v5.1.0/bio/gatk/haplotypecaller"


#
# rule tabix_all_variants:
#     input:
#         "results/called/{sample}.{contig}.g.vcf.gz",
#     output:
#         "results/called/{sample}.{contig}.g.vcf.gz.tbi",
#     log:
#         "logs/called/tabix/{sample}.{contig}.log",
#     params:
#         "-p vcf",
#     wrapper:
#         "v5.1.0/bio/tabix/index"


rule combine_calls:
    input:
        ref="resources/chromosomes.fasta",
        tabix=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz.tbi", sample=samples.index
        ),
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf=temp("results/called/all.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    threads: 2
    params:
        java_opts="-XX:ParallelGCThreads=2",
    resources:
        mem_mb=64000,
    wrapper:
        "v5.1.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/chromosomes.fasta",
        gvcf="results/called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    resources:
        mem_mb=64000,
    wrapper:
        "v5.1.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    threads: 4
    params:
        java_opts="-XX:ParallelGCThreads=4",
    resources:
        mem_mb=64000,
    wrapper:
        "v5.1.0/bio/picard/mergevcfs"
