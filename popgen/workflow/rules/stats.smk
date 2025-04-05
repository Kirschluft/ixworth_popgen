rule worst_consequence_stats:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/consequence_worst.tsv.gz",
    log:
        "logs/consequence_worst.log",
    conda:
        "../envs/stats.yaml"
    shell:
        """
        bcftools +split-vep -d -s worst \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Consequence\t%IMPACT\t%Gene\t%Amino_acids\t%SIFT\n' \
        {input} -Oz -o {output} > {log} 2>&1 
        """


rule all_consequence_stats:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/consequence_all.tsv.gz",
    log:
        "logs/stats/consequence_all.log",
    conda:
        "../envs/stats.yaml"
    shell:
        """
        bcftools +split-vep -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%Consequence\t%IMPACT\t%Gene\t%Amino_acids\t%SIFT\n' \
        {input} -Oz -o {output} > {log} 2>&1 
        """


rule vcf_stats:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/vcf_stats.tsv.gz",
    log:
        "logs/stats/vcf_stats.log",
    conda:
        "../envs/stats.yaml"
    shell:
        """
        bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/AF\t%INFO/DP\t%INFO/BaseQRankSum\t%INFO/ExcessHet\t%INFO/FS\t%INFO/InbreedingCoeff\t%INFO/MLEAC\t%INFO/MLEAF\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/QD\t%INFO/ReadPosRankSum\t%INFO/SOR\t%INFO/AN\n' \
        {input} | gzip -c > {output} 2> {log}
        """


rule intersect_populations:
    input:
        vcf=expand(
            "results/vcfs/populations/{population}.vcf.gz",
            population=list(POPULATIONS.keys()),
        ),
        index=expand(
            "results/vcfs/populations/{population}.vcf.gz.csi",
            population=list(POPULATIONS.keys()),
        ),
    output:
        "results/stats/intersection.txt",
    log:
        "logs/intersect.log",
    conda:
        "../envs/stats.yaml"
    shell:
        "vcf-compare {input.vcf} > {output} 2> {log}"


rule visualise_intersection:
    input:
        "results/stats/intersection.txt",
    output:
        "results/stats/intersection_upset.png",
    log:
        "logs/stats/intersection_upset.log",
    conda:
        "../envs/stats.yaml"
    shell:
        """
        Rscript workflow/scripts/intersect.R --out {output} --file {input} \
        > {log} 2>&1
        """


rule vcftools_het:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/stats/het/{population}.het",
    log:
        "logs/stats/het_{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --het --out results/stats/het/{wildcards.population} > {log} 2>&1"


rule vcftools_het_all:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.het",
    log:
        "logs/stats/het_populations.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --het --out results/stats/populations > {log} 2>&1"


rule vcftools_depth:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.idepth",
    log:
        "logs/stats/depth.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --depth --out results/stats/populations > {log} 2>&1"


rule vcftools_site_mean_depth:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.ldepth.mean",
    log:
        "logs/stats/depth_site_mean.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --site-mean-depth --out results/stats/populations > {log} 2>&1"


rule vcftools_site_mean_quality:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.lqual",
    log:
        "logs/stats/site_quality.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --site-quality --out results/stats/populations > {log} 2>&1"


rule vcftools_relatedness:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.relatedness",
    log:
        "logs/stats/relatedness.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --relatedness --out results/stats/populations > {log} 2>&1"


rule vcftools_tstv:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/stats/tstv/{population}.TsTv.summary",
    log:
        "logs/stats/tstv_{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --TsTv-summary --out results/stats/tstv/{wildcards.population} > {log} 2>&1"


rule vcftools_tstv_all:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.TsTv.summary",
    log:
        "logs/stats/tstv_populations.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --TsTv-summary --out results/stats/populations > {log} 2>&1"


rule vcftools_af:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/stats/freq/{population}.frq",
    log:
        "logs/stats/freq_{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --freq --out results/stats/freq/{wildcards.population} > {log} 2>&1"


rule vcftools_af_all:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.frq",
    log:
        "logs/stats/tstv.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --freq --out results/stats/populations > {log} 2>&1"


rule vcftools_missingness:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/stats/lmiss/{population}.lmiss",
    log:
        "logs/stats/lmiss_{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-site --out results/stats/lmiss/{wildcards.population} > {log} 2>&1"


rule vcftools_missingness_all:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/stats/populations.lmiss",
    log:
        "logs/stats/lmiss_populations.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-site --out results/stats/populations > {log} 2>&1"


rule het_stats:
    input:
        expand(
            "results/stats/het/{population}.het", population=list(POPULATIONS.keys())
        ),
    output:
        "results/stats/het/het_F_violinplot.png",
    log:
        "logs/stats/het_plots.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "Rscript workflow/scripts/het.R --out results/stats/het/het --files {input} > {log} 2>&1"
