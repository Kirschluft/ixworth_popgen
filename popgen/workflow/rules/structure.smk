rule convert_plink:
    input:
        "results/vcfs/annotated.vcf.gz",
    output:
        bed="results/plink/filtered.bed",
        bim="results/plink/filtered.bim",
        fam="results/plink/filtered.fam",
    log:
        "logs/plink/convert_plink.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --vcf {input} --make-bed --chr-set {N_CHROMOSOMES} \
            --allow-extra-chr --out results/plink/filtered > {log} 2>&1
        """


rule add_population_plink:
    input:
        bed="results/plink/filtered.bed",
        bim="results/plink/filtered.bim",
        fam="results/plink/filtered.fam",
    output:
        ped="results/plink/filtered_population.ped",
        m="results/plink/filtered_population.map",
    log:
        "logs/plink/add_population_plink.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --bfile results/plink/filtered \
         --chr-set {N_CHROMOSOMES} --update-ids {config[populations_map]} \
         --allow-extra-chr --recode --out results/plink/filtered_population > {log} 2>&1
        """


rule prune:
    input:
        ped="results/plink/filtered_population.ped",
        m="results/plink/filtered_population.map",
    output:
        pruned="results/plink/pruned.prune.in",
    log:
        "logs/plink/prune.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --file results/plink/filtered_population --indep-pairwise {config[prune][pca_admixture]} \
            --chr-set {N_CHROMOSOMES} --allow-extra-chr \
            --out results/plink/pruned > {log} 2>&1
        """


rule pruned_vcf:
    input:
        pruned="results/plink/pruned.prune.in",
    output:
        "results/vcfs/pruned.vcf",
    log:
        "logs/plink/pruned_vcf.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --file results/plink/filtered_population --extract {input.pruned} \
            --recode vcf --chr-set {N_CHROMOSOMES} --allow-extra-chr \
            --out results/vcfs/pruned > {log} 2>&1
        """


rule pruned_bed:
    input:
        ped="results/plink/filtered_population.ped",
        m="results/plink/filtered_population.map",
        pruned="results/plink/pruned.prune.in",
    output:
        bed="results/plink/pruned.bed",
        bim="results/plink/pruned.bim",
        fam="results/plink/pruned.fam",
    log:
        "logs/plink/thin.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        plink --file results/plink/filtered_population --extract {input.pruned} \
            --chr-set {N_CHROMOSOMES} --allow-extra-chr --make-bed \
            --out results/plink/pruned > {log} 2>&1
        """


rule convert_vcf_gds:
    input:
        "results/vcfs/pruned.vcf",
    output:
        "results/vcfs/pruned.gds",
    log:
        "logs/pca/convert.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        Rscript workflow/scripts/vcf2gds.R --vcf {input} \
            --out results/vcfs/pruned.gds > {log} 2>&1
         """


rule pca:
    input:
        "results/vcfs/pruned.gds",
    output:
        pca="results/pca/pca.png",
        pca5="results/pca/pca_top5.png",
    log:
        "logs/pca/pca.log",
    conda:
        "../envs/structure.yaml"
    shell:
        """
        Rscript workflow/scripts/pca.R --file {input} \
            --out results/pca/ > {log} 2>&1
         """


rule roh:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/roh/{population}_regions.txt",
    log:
        "results/roh/{population}_regions.log",
    conda:
        "../envs/structure.yaml"
    threads: 16
    shell:
        """
        bcftools roh --AF-dflt 0.4 --estimate-AF - --threads {threads} {input} > {output} 2> {log}
        """


rule roh_filter:
    input:
        "results/roh/{population}_regions.txt",
    output:
        "results/roh/{population}.txt",
    log:
        "results/roh/{population}.log",
    shell:
        """
        grep '^RG' {input} > {output} 2> {log}
        """


rule roh_stats:
    input:
        expand(
            "results/roh/{population}.txt",
            population=list(POPULATIONS.keys()),
        ),
    output:
        roh="results/roh/roh_count_mean.png",
    log:
        "logs/roh/roh_stats.log",
    conda:
        "../envs/structure.yaml"
    params:
        length=config["genome_length"],
    shell:
        """
        Rscript workflow/scripts/roh.R --input {input} --length {params.length} --out results/roh > {log} 2>&1
        """


rule admixture_k:
    input:
        bed="results/plink/pruned.bed",
        bim="results/plink/pruned.bim",
        fam="results/plink/pruned.fam",
    output:
        P="results/admixture/pruned.{k}.P",
        Q="results/admixture/pruned.{k}.Q",
    log:
        "logs/admixture/k{k}.log",
    conda:
        "../envs/structure.yaml"
    threads: workflow.cores
    params:
        cv=lambda wc: config["admixture"]["c"],
        seed=lambda wc: int(wc.k) * 42,
    shell:
        """
        mkdir -p results/admixture logs/admixture
        admixture --cv={params.cv} -j{threads} {input.bed} {wildcards.k} --seed={params.seed} > {log} 2>&1
        mv pruned.{wildcards.k}.P pruned.{wildcards.k}.Q results/admixture/
        """


rule collect_admixture_cv:
    input:
        logs=expand(
            "logs/admixture/k{k}.log",
            k=range(config["admixture"]["k"], config["admixture"]["K"] + 1),
        ),
    output:
        "results/admixture/cv.txt",
    shell:
        "grep -h 'CV' {input.logs} > {output}"


rule visualize_admixture:
    input:
        cv="results/admixture/cv.txt",
        fam="results/plink/pruned.fam",
        runs=expand(
            "results/admixture/pruned.{k}.P",
            k=range(config["admixture"]["k"], config["admixture"]["K"] + 1),
        ),
    output:
        crossval="results/admixture/cross_validation.png",
        admix="results/admixture/admixture.png",
    log:
        "logs/admixture/plot.log",
    conda:
        "../envs/structure.yaml"
    shell:
        r"""
        Rscript workflow/scripts/admixture.R --prefix results/admixture/pruned \
            --fam {input.fam} --cv {input.cv} --out results/admixture/ > {log} 2>&1
        """
