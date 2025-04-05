rule vcftools_fst_windows:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        population1=lambda wildcards: POPULATIONS[wildcards.pop1],
        population2=lambda wildcards: POPULATIONS[wildcards.pop2],
    output:
        "results/fst/windows/{pop1}_vs_{pop2}.windowed.weir.fst",
    log:
        "logs/fst/windows/{pop1}_vs_{pop2}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --weir-fst-pop {input.population1} \
                 --weir-fst-pop {input.population2} \
                 --fst-window-size {config[windows][size]} \
                 --fst-window-step {config[windows][step]} \
                 --out results/fst/windows/{wildcards.pop1}_vs_{wildcards.pop2} \
                 > {log} 2>&1
        """


rule vcftools_fst_sites:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        population1=lambda wildcards: POPULATIONS[wildcards.pop1],
        population2=lambda wildcards: POPULATIONS[wildcards.pop2],
    output:
        "results/fst/sites/{pop1}_vs_{pop2}.weir.fst",
    log:
        "logs/fst/sites/{pop1}_vs_{pop2}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --weir-fst-pop {input.population1} \
                 --weir-fst-pop {input.population2} \
                 --out results/fst/sites/{wildcards.pop1}_vs_{wildcards.pop2} \
                 > {log} 2>&1
        """


rule fst_stats:
    input:
        populations=expand(
            "results/fst/windows/{pop1}_vs_{pop2}.windowed.weir.fst",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
        genes="results/resources/genes.gff3.gz",
    output:
        sum="results/fst/summary.csv",
        box="results/fst/boxplot.png",
        violin="results/fst/violinplot.png",
        sum_chr="results/fst/chr_summary.csv",
        bar_chr="results/fst/bar_chr.png",
        dense_chr="results/fst/dense_chr.png",
        heatmap="results/fst/heatmap.png",
        heatmap_chr="results/fst/heatmap_chr.png",
        sweeps=expand(
            "results/fst/{pop1}_vs_{pop2}_sweeps.csv",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
        freqs=expand(
            "results/fst/{pop1}_vs_{pop2}_sweeps_freq.png",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
        hist_snp=expand(
            "results/fst/{pop1}_vs_{pop2}_nsnps.png",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
        hist_fst=expand(
            "results/fst/{pop1}_vs_{pop2}_hist.png",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
        manhattan=expand(
            "results/fst/{pop1}_vs_{pop2}_sweeps.png",
            zip,
            pop1=FST_CONTRAST1,
            pop2=FST_CONTRAST2,
        ),
    log:
        "logs/fst/stats.log",
    conda:
        "../envs/vcftools.yaml"
    threads: 16
    shell:
        """Rscript workflow/scripts/fst.R --genes {input.genes} \
            --files {input.populations} \
            --cores {threads} \
            --chromosomes {CHROMOSOMES} \
            --out results/fst/ > {log} 2>&1
         """


rule vcftools_pi:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        population=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/pi/{population}.windowed.pi",
    log:
        "logs/pi/{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --keep {input.population} \
                 --window-pi {config[windows][size]} \
                 --window-pi-step {config[windows][step]} \
                 --out results/pi/{wildcards.population} \
                 > {log} 2>&1
        """


rule pi_stats:
    input:
        expand(
            "results/pi/{population}.windowed.pi", population=list(POPULATIONS.keys())
        ),
    output:
        sum="results/pi/summary.csv",
        box="results/pi/boxplot.png",
        violin="results/pi/violinplot.png",
        sum_chr="results/pi/chr_summary.csv",
        bar_chr="results/pi/bar_chr.png",
        dense_chr="results/pi/dense_chr.png",
        heatmap_chr="results/pi/heatmap_chr.png",
        hist_pi=expand(
            "results/pi/{population}_hist.png", population=list(POPULATIONS.keys())
        ),
        manhattan=expand(
            "results/pi/{population}_manhattan.png",
            population=list(POPULATIONS.keys()),
        ),
    log:
        "logs/pi/stats.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "Rscript workflow/scripts/pi.R --out results/pi/ --files {input} > {log} 2>&1"


rule vcftools_tajimasd:
    input:
        vcf="results/vcfs/filtered.vcf.gz",
        population=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/tajimas/{population}.Tajima.D",
    log:
        "logs/tajimas/{population}.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --keep {input.population} \
                 --out results/tajimas/{wildcards.population} \
                 --TajimaD {config[windows][size]} \
                 > {log} 2>&1
        """


rule tajimasd_stats:
    input:
        expand(
            "results/tajimas/{population}.Tajima.D",
            population=list(POPULATIONS.keys()),
        ),
    output:
        sum="results/tajimas/summary.csv",
        box="results/tajimas/boxplot.png",
        violin="results/tajimas/violinplot.png",
        sum_chr="results/tajimas/chr_summary.csv",
        bar_chr="results/tajimas/bar_chr.png",
        dense_chr="results/tajimas/dense_chr.png",
        heatmap_chr="results/tajimas/heatmap_chr.png",
        hist_pi=expand(
            "results/tajimas/{population}_hist.png",
            population=list(POPULATIONS.keys()),
        ),
        manhattan=expand(
            "results/tajimas/{population}_manhattan.png",
            population=list(POPULATIONS.keys()),
        ),
    log:
        "logs/tajimas/stats.log",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "Rscript workflow/scripts/tajima.R --out results/tajimas/ --files {input} > {log} 2>&1"
