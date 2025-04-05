rule index_result_vcf:
    input:
        VCF_FILE,
    output:
        f"{VCF_FILE}.csi",
    log:
        "logs/bcftools_result_index.log",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule filter_result_snps:
    input:
        vcf=VCF_FILE,
        index=f"{VCF_FILE}.csi",
    output:
        "results/vcfs/snps.vcf.gz",
    log:
        "logs/bcftools_result_snps.log",
    threads: 16
    params:
        extra=f"-f PASS -r {','.join(CHROMOSOMES)} --max-alleles 2 --exclude-types indels",
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule index_snps_vcf:
    input:
        "results/vcfs/snps.vcf.gz",
    output:
        "results/vcfs/snps.vcf.gz.csi",
    log:
        "logs/bcftools_snps_index.log",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule filter_result_missing:
    input:
        vcf="results/vcfs/snps.vcf.gz",
        index="results/vcfs/snps.vcf.gz.csi",
    output:
        "results/vcfs/filtered.vcf.gz",
    log:
        "logs/bcftools_result_filtered.log",
    params:
        extra="-i 'F_MISS <= 0.5 && MAC >= 2'",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule index_filtered_result:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/vcfs/filtered.vcf.gz.csi",
    log:
        "logs/bcftools_filtered_result_index.log",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule annotate_vcf:
    input:
        "results/vcfs/filtered.vcf.gz",
    output:
        "results/vcfs/annotated.vcf.gz",
    log:
        "logs/bcftools/annotate_vcf.log",
    conda:
        "../envs/structure.yaml"
    threads: 16
    shell:
        """
        bcftools annotate --set-id 'VAR\_%CHROM\_%POS' \
            -Oz -o {output} {input} --threads {threads} > {log} 2>&1
        """


rule index_annotated_result:
    input:
        "results/vcfs/annotated.vcf.gz",
    output:
        "results/vcfs/annotated.vcf.gz.csi",
    log:
        "logs/bcftools_annotated_result_index.log",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"


rule split_population:
    input:
        "results/vcfs/annotated.vcf.gz",
        index="results/vcfs/annotated.vcf.gz.csi",
        samples=lambda wildcards: POPULATIONS[wildcards.population],
    output:
        "results/vcfs/populations/{population}.vcf.gz",
    log:
        "logs/split_population/{population}.log",
    threads: 16
    params:
        extra="-c 1",
    wrapper:
        "v5.1.0/bio/bcftools/view"


rule index_population:
    input:
        "results/vcfs/populations/{population}.vcf.gz",
    output:
        "results/vcfs/populations/{population}.vcf.gz.csi",
    log:
        "logs/index_population/{population}.log",
    threads: 16
    wrapper:
        "v5.1.0/bio/bcftools/index"
