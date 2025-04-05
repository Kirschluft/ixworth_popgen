rule bcftools_stats:
    input:
        "results/filtered/all.vcf.gz",
    output:
        "results/stats/bcftools_stats.txt",
    log:
        "logs/stats/bcftools.log",
    wrapper:
        "v5.1.0/bio/bcftools/stats"


rule bcftools_plot:
    input:
        "results/stats/bcftools_stats.txt",
    output:
        directory("results/stats/bcftools_plots/"),
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/stats/bcftools.log",
    shell:
        """
        plot-vcfstats -p {output} {input} > {log} 2>&1
        """


# Workaround ...
rule collect_metrics:
    input:
        filtered="results/filtered/all.vcf.gz",
        known="resources/variation.chromosomes.reheadered.vcf.gz",
        dict="resources/chromosomes.dict",
    output:
        metrics=directory("results/stats/picard_metrics"),
    conda:
        "../envs/picard.yaml"
    log:
        "logs/picard/metric_collection/all.metrics.log",
    shell:
        """
        picard CollectVariantCallingMetrics \
            -I {input.filtered} \
            --DBSNP {input.known} \
            -SD {input.dict} \
            -O {output.metrics} > {log} 2>&1
        mkdir -p {output.metrics}
        mv results/stats/picard_metrics.* {output.metrics}
        """
