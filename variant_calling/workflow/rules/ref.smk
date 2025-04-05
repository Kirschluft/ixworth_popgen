rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v5.1.0/bio/reference/ensembl-sequence"


rule get_chromosomes:
    input:
        "resources/genome.fasta",
    output:
        "resources/chromosomes.fasta",
    log:
        "logs/get-chromosomes.log",
    params:
        command="grep",
        extra="-p '^[1-9WZXY]' -r",
    threads: 2
    cache: True
    wrapper:
        "v5.1.0/bio/seqkit"


checkpoint genome_faidx:
    input:
        "resources/chromosomes.fasta",
    output:
        "resources/chromosomes.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v5.1.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/chromosomes.fasta",
    output:
        "resources/chromosomes.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/chromosomes.fasta.fai",
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "v5.1.0/bio/reference/ensembl-variation"


rule tabix_known_variants:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "v5.1.0/bio/tabix/index"


# workaround for picard
rule filter_known_variation:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.chromosomes.vcf.gz",
    log:
        "logs/variation.filtered.log",
    params:
        # will only work for common chr names
        extra="-i \"CHROM ~ '^chr[1-9][0-9]*$' || CHROM ~ '^[1-9][0-9]*$' || CHROM == 'X' || CHROM == 'Y' || CHROM == 'W' || CHROM == 'Z'\"",
    cache: True
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule reheader_known_variation:
    input:
        vcf="resources/variation.chromosomes.vcf.gz",
        fai="resources/chromosomes.fasta.fai",
    output:
        "resources/variation.chromosomes.reheadered.vcf.gz",
    log:
        "logs/variation.filtered.reheader.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        """
      bcftools reheader -f {input.fai} {input.vcf} -o {output} &> {log}
      """


rule tabix_filtered_variants:
    input:
        "resources/variation.chromosomes.vcf.gz",
    output:
        "resources/variation.chromosomes.vcf.gz.tbi",
    log:
        "logs/tabix/variation.filtered.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "v5.1.0/bio/tabix/index"


rule bwa_index:
    input:
        "resources/chromosomes.fasta",
    output:
        multiext("resources/chromosomes", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v5.1.0/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "v5.1.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "v5.1.0/bio/vep/plugins"
