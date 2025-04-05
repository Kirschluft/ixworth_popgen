rule fix_ploidy:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/vcfs/filtered_diploid.vcf.gz",
    log:
        "logs/vcfs/fix_ploidy.log",
    conda:
        "../envs/haplo.yaml"
    shell:
        "bcftools +fixploidy {input} -- -f 2 | bgzip > {output} 2> {log}"


rule index_fixed_vcf:
    input:
        "results/vcfs/filtered_diploid.vcf.gz",
    output:
        "results/vcfs/filtered_diploid.vcf.gz.csi",
    log:
        "logs/vcfs/index_fixed_vcf.log",
    conda:
        "../envs/haplo.yaml"
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule beagle_phasing:
    input:
        vcf="results/vcfs/filtered_diploid.vcf.gz",
    output:
        "results/phased/filtered_phased.vcf.gz",
    log:
        "logs/phased/beagle_phasing.log",
    conda:
        "../envs/haplo.yaml"
    threads: 64
    params:
        mem=config["beagle"]["mem"],
        iterations=config["beagle"]["iterations"],
        burnin=config["beagle"]["burnin"],
        extra=config["beagle"]["extra"],
    shell:
        """
        beagle -Xmx{params.mem} gt={input} seed=42 \
        out=results/phased/filtered_phased nthreads={threads} \
        iterations={params.iterations} burnin={params.burnin} {params.extra} > {log} 2>&1
        """


rule index_phased_vcf:
    input:
        "results/phased/filtered_phased.vcf.gz",
    output:
        "results/phased/filtered_phased.vcf.gz.csi",
    log:
        "logs/phased/index_phased_vcf.log",
    conda:
        "../envs/haplo.yaml"
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule get_ancestral_allele:
    input:
        vcf="results/phased/filtered_phased.vcf.gz",
        index="results/phased/filtered_phased.vcf.gz.csi",
    output:
        "results/phased/ancestral_alleles.txt.gz",
    log:
        "logs/phased/get_ancestral_allele.log",
    conda:
        "../envs/haplo.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%REF\n' {input.vcf} | bgzip > {output} 2> {log}"


rule index_ancestral_allele:
    input:
        "results/phased/ancestral_alleles.txt.gz",
    output:
        "results/phased/ancestral_alleles.txt.gz.tbi",
    log:
        "logs/phased/index_ancestral_allele.log",
    conda:
        "../envs/haplo.yaml"
    params:
        extra="-s1 -b2 -e2 -h",
    wrapper:
        "v5.1.0/bio/tabix/index"


rule annotate_ancestral_allele:
    input:
        vcf="results/phased/filtered_phased.vcf.gz",
        vcf_index="results/phased/filtered_phased.vcf.gz.csi",
        ancestral="results/phased/ancestral_alleles.txt.gz",
        ancestral_index="results/phased/ancestral_alleles.txt.gz.tbi",
    output:
        vcf="results/phased/final_phased.vcf.gz",
        header="results/phased/ancestral_alleles.header",
    log:
        "logs/phased/annotate_ancestral_allele.log",
    conda:
        "../envs/haplo.yaml"
    shell:
        """
        echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">' > {output.header}
        bcftools annotate -a {input.ancestral} -h {output.header} \
            -c CHROM,POS,INFO/AA {input.vcf} -Oz -o {output.vcf}
        """


rule index_ancestral_vcf:
    input:
        "results/phased/final_phased.vcf.gz",
    output:
        "results/phased/final_phased.vcf.gz.csi",
    log:
        "logs/phased/index_ancestral_vcf.log",
    conda:
        "../envs/haplo.yaml"
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule split_phased_vcf:
    input:
        "results/phased/final_phased.vcf.gz",
        index="results/phased/final_phased.vcf.gz.csi",
        samples=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/phased/{population}.phased.vcf.gz",
    log:
        "logs/phased/{population}.log",
    conda:
        "../envs/haplo.yaml"
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule ihs:
    input:
        phased="results/phased/{population}.phased.vcf.gz",
    output:
        ihs="results/ihs/{population}_ihs.csv.gz",
        freq="results/ihs/{population}_freq.csv.gz",
    log:
        "logs/ihs/runs/{population}.log",
    container:
        "workflow/envs/rehh.sif"
    threads: 16
    shell:
        """
        Rscript workflow/scripts/ihs.R --file {input.phased} --chromosomes {CHROMOSOMES} --cores {threads} \
            --out results/ihs/{wildcards.population}  > {log} 2>&1
        """


rule ihs_stats:
    input:
        ihs="results/ihs/{population}_ihs.csv.gz",
        genes="results/resources/genes.gff3.gz",
    output:
        hist="results/ihs/{population}_hist.png",
        regions="results/ihs/{population}_regions.csv",
        sweeps="results/ihs/{population}_sweeps.csv",
        freqs="results/ihs/{population}_sweeps_freq.png",
        manhattan="results/ihs/{population}_sweeps.png",
    log:
        "logs/ihs/stats/{population}.log",
    container:
        "workflow/envs/rehh.sif"
    shell:
        """
        Rscript workflow/scripts/ihs_stats.R --input {input.ihs} --genes {input.genes} \
            --window_size {config[windows][size]} --window_step {config[windows][step]} \
            --chromosomes {CHROMOSOMES} --out results/ihs/{wildcards.population}  > {log} 2>&1
        """


rule xp_ehh:
    input:
        phased1="results/phased/{pop1}.phased.vcf.gz",
        phased2="results/phased/{pop2}.phased.vcf.gz",
    output:
        "results/xp_ehh/{pop1}_vs_{pop2}_xp_ehh.csv.gz",
    log:
        "logs/xp_ehh/runs/{pop1}_vs_{pop2}.log",
    container:
        "workflow/envs/rehh.sif"
    threads: 16
    shell:
        """
        Rscript workflow/scripts/xp_ehh.R --file1 {input.phased1} --file2 {input.phased2} \
            --chromosomes {CHROMOSOMES} --cores {threads} \
            --out results/xp_ehh/{wildcards.pop1}_vs_{wildcards.pop2} > {log} 2>&1
        """


rule xp_ehh_stats:
    input:
        xp_ehh="results/xp_ehh/{pop1}_vs_{pop2}_xp_ehh.csv.gz",
        genes="results/resources/genes.gff3.gz",
    output:
        hist="results/xp_ehh/{pop1}_vs_{pop2}_hist.png",
        regions="results/xp_ehh/{pop1}_vs_{pop2}_regions.csv",
        sweeps="results/xp_ehh/{pop1}_vs_{pop2}_sweeps.csv",
        freqs="results/xp_ehh/{pop1}_vs_{pop2}_sweeps_freq.png",
        manhattan="results/xp_ehh/{pop1}_vs_{pop2}_sweeps.png",
    log:
        "logs/xp_ehh/stats/{pop1}_vs_{pop2}.log",
    container:
        "workflow/envs/rehh.sif"
    shell:
        """
        Rscript workflow/scripts/xp_ehh_stats.R --input {input.xp_ehh} --genes {input.genes} \
            --window_size {config[windows][size]} --window_step {config[windows][step]}  \
            --chromosomes {CHROMOSOMES} --out results/xp_ehh/{wildcards.pop1}_vs_{wildcards.pop2} > {log} 2>&1
        """


rule fst_go:
    input:
        fst_sweeps="results/fst/{pop1}_vs_{pop2}_sweeps.csv",
        qtls=config["qtls"],
    output:
        "results/function/{pop1}_vs_{pop2}_fst_go_terms.csv",
    log:
        "logs/function/{pop1}_vs_{pop2}_fst.log",
    conda:
        "../envs/go.yaml"
    shell:
        """
        Rscript workflow/scripts/go_function.R --file {input.fst_sweeps} \
            --gprofiler_version https://biit.cs.ut.ee/gprofiler_archive3/e106_eg53_p16/ --organism ggallus --qtls {input.qtls} \
            --out results/function/{wildcards.pop1}_vs_{wildcards.pop2}_fst > {log} 2>&1
        """


rule xp_ehh_go:
    input:
        xp_ehh_sweeps="results/xp_ehh/{pop1}_vs_{pop2}_sweeps.csv",
        qtls=config["qtls"],
    output:
        "results/function/{pop1}_vs_{pop2}_xp_ehh_go_terms.csv",
    log:
        "logs/function/{pop1}_vs_{pop2}_xp_ehh.log",
    conda:
        "../envs/go.yaml"
    shell:
        """
        Rscript workflow/scripts/go_function.R --file {input.xp_ehh_sweeps} \
            --gprofiler_version https://biit.cs.ut.ee/gprofiler_archive3/e106_eg53_p16/ --organism ggallus --qtls {input.qtls} \
            --out results/function/{wildcards.pop1}_vs_{wildcards.pop2}_xp_ehh > {log} 2>&1
        """


rule ihs_go:
    input:
        ihs_sweeps="results/ihs/{population}_sweeps.csv",
        qtls=config["qtls"],
    output:
        "results/function/{population}_ihs_go_terms.csv",
    log:
        "logs/function/{population}_ihs.log",
    conda:
        "../envs/go.yaml"
    shell:
        """
        Rscript workflow/scripts/go_function.R --file {input.ihs_sweeps} \
            --gprofiler_version https://biit.cs.ut.ee/gprofiler_archive3/e106_eg53_p16/ --organism ggallus --qtls {input.qtls} \
            --out results/function/{wildcards.population}_ihs > {log} 2>&1
        """
