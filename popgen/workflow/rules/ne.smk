rule convert_population_plink:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        bed="results/plink/populations/{population}.bed",
        bim="results/plink/populations/{population}.bim",
        fam="results/plink/populations/{population}.fam",
    log:
        "logs/plink/{population}/convert_plink.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --chr-set {N_CHROMOSOMES} \
            --allow-extra-chr --out results/plink/populations/{wildcards.population} > {log} 2>&1
        """


rule prune_population_snps:
    input:
        bed="results/plink/populations/{population}.bed",
        bim="results/plink/populations/{population}.bim",
        fam="results/plink/populations/{population}.fam",
    output:
        "results/plink/populations/{population}.prune.in",
    log:
        "logs/ne/pruning/{population}.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --bfile results/plink/populations/{wildcards.population} --indep-pairwise {config[prune][ne][ld]} \
            --chr-set {N_CHROMOSOMES} --allow-extra-chr \
            --out results/plink/populations/{wildcards.population} > {log} 2>&1
        """


rule convert_pruned_population_snps:
    input:
        "results/plink/populations/{population}.prune.in",
    output:
        "results/ne/{population}_pruned.vcf",
    log:
        "logs/ne/converting/{population}.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --bfile results/plink/populations/{wildcards.population} --extract {input} \
            --recode vcf --chr-set {N_CHROMOSOMES} --bp-space {config[prune][ne][bp_space]} \
            --maf {config[prune][ne][maf]} \
            --allow-extra-chr --out results/ne/{wildcards.population}_pruned > {log} 2>&1
        """


rule index_pruned_population:
    input:
        "results/ne/{population}_pruned.vcf",
    output:
        "results/ne/{population}_pruned.vcf.csi",
    log:
        "logs/index_population/pruned/{population}.log",
    threads: 16
    conda:
        "../envs/structure.yaml"
    wrapper:
        "v5.1.0/bio/bcftools/index"


# Depending on dartR, this might change in the future
rule get_neestimator:
    output:
        "results/resources/Neestimator/neestimator_linux.zip",
        "results/resources/Neestimator/Ne2-1L",
    container:
        "workflow/envs/ne.sif"
    shell:
        """
        wget -q -N https://github.com/green-striped-gecko/dartRverse/raw/main/binaries/neestimator_linux.zip -O {output[0]}
        unzip {output[0]} -d results/resources/
        chmod +x {output[1]}
        """


rule run_ne_estimator:
    input:
        nee="results/resources/Neestimator/Ne2-1L",
        vcf="results/ne/{population}_pruned.vcf",
    output:
        "results/ne/ne_estimator/{population}_contemporary_ne.csv",
    log:
        "logs/ne_estimator/{population}.log",
    container:
        "workflow/envs/ne.sif"
    shell:
        """
    Rscript workflow/scripts/ne.R --file {input.vcf} --path {input.nee} \
    --nchr {N_CHROMOSOMES} --out {output} > {log} 2>&1
    """


rule split_smc:
    input:
        vcf="results/vcfs/populations/{population}.vcf.gz",
        index="results/vcfs/populations/{population}.vcf.gz.csi",
    output:
        "results/ne/smcpp/{population}/chr{chr}.smc.gz",
    log:
        "logs/smcpp/{population}/{chr}.log",
    conda:
        "../envs/smcpp.yaml"
    params:
        indiv=lambda wildcards: ",".join(
            open(POPULATIONS[wildcards.population]).read().splitlines()
        ),
    shell:
        """
    smc++ vcf2smc {input.vcf} results/ne/smcpp/{wildcards.population}/chr{wildcards.chr}.smc.gz \
    {wildcards.chr} {wildcards.population}:{params.indiv} > {log} 2>&1
    """


rule run_smcpp:
    input:
        vcfs=expand(
            "results/ne/smcpp/{population}/chr{chr}.smc.gz",
            population="{population}",
            chr=CHROMOSOMES,
        ),
    output:
        "results/ne/smcpp/{population}/estimates/model.final.json",
    log:
        "logs/smcpp/estimates/{population}.log",
    conda:
        "../envs/smcpp.yaml"
    threads: 16
    shell:
        """
    smc++ cv --folds {config[ne][folds]} --cores {threads} {config[ne][mu]} \
    -o results/ne/smcpp/{wildcards.population}/estimates/ results/ne/smcpp/{wildcards.population}/chr*.smc.gz > {log} 2>&1
    """


rule gather_smcpp:
    input:
        pops=expand(
            "results/ne/smcpp/{population}/estimates/model.final.json",
            population=list(POPULATIONS.keys()),
        ),
    output:
        csv="results/ne/smcpp/historical_ne.csv",
        png="results/ne/smcpp/historical_ne.png",
    log:
        "logs/smcpp/plot.log",
    conda:
        "../envs/smcpp.yaml"
    shell:
        """
    smc++ plot -g {config[ne][generation_time]} -c {output.png} {input.pops} > {log} 2>&1
    """
