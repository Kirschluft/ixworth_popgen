rule split_vcf:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        index="results/vcfs/filtered.vcf.gz.csi",
    output:
        "results/vcfs/chr/chr{chr}.vcf.gz",
    log:
        "logs/ld/chr{chr}_split.log",
    params:
        extra="-r {chr}",
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule run_ld_decay_pop:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        pop_file=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/ld/pop/{population}.stat.gz",
    log:
        "logs/ld/pop/{population}.log",
    conda:
        "../envs/ld.yaml"
    shell:
        """
        PopLDdecay -InVCF {input.vcf} -OutStat results/ld/pop/{wildcards.population} \
          -SubPop {input.pop_file} -MaxDist 500 -Het 1.0 -Miss 0 -MAF 0 -OutType 2 > {log} 2>&1
        """


rule run_ld_decay_chr:
    input:
        vcf="results/vcfs/chr/chr{chr}.vcf.gz",
        pop_file=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/ld/chr/{population}_chr{chr}.stat.gz",
    log:
        "logs/ld/chr/{population}_chr{chr}.log",
    conda:
        "../envs/ld.yaml"
    shell:
        """
        PopLDdecay -InVCF {input.vcf} -OutStat results/ld/chr/{wildcards.population}_chr{wildcards.chr} -SubPop {input.pop_file} \
        -MaxDist 500 -Het 1.0 -Miss 0 -MAF 0 -OutType 2 > {log} 2>&1 
        """


rule visualize_ld_decay_pop:
    input:
        expand(
            "results/ld/pop/{population}.stat.gz", population=list(POPULATIONS.keys())
        ),
    output:
        out1="results/ld/populations.png",
        out2="results/ld/populations_close.png",
    log:
        "logs/ld/visualisation_populations.log",
    conda:
        "../envs/ld.yaml"
    shell:
        "Rscript workflow/scripts/ld.R --out results/ld/populations --files {input} > {log} 2>&1"


rule visualize_ld_decay_chr:
    input:
        lambda wildcards: expand(
            "results/ld/chr/{population}_chr{chr}.stat.gz",
            population=wildcards.population,
            chr=CHROMOSOMES,
        ),
    output:
        out1="results/ld/{population}_chr.png",
        out2="results/ld/{population}_chr_close.png",
    log:
        "logs/ld/visualisation_chr_{population}.log",
    conda:
        "../envs/ld.yaml"
    shell:
        "Rscript workflow/scripts/ld.R --out results/ld/{wildcards.population}_chr --files {input} --column Chromosome > {log} 2>&1"
