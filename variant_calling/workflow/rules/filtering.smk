rule select_calls:
    input:
        ref="resources/chromosomes.fasta",
        vcf="results/genotyped/all.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.vcf.gz"),
    params:
        extra=get_vartype_arg,
    log:
        "logs/gatk/selectvariants/{vartype}.log",
    wrapper:
        "v5.1.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/chromosomes.fasta",
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.hardfiltered.vcf.gz"),
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    wrapper:
        "v5.1.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="results/filtered/all.{vartype}.vcf.gz",
        ref="resources/chromosomes.fasta",
        dict="resources/chromosomes.dict",
        ensembl="resources/variation.noiupac.vcf.gz",
        ensembl_idx="resources/variation.noiupac.vcf.gz.tbi",
    output:
        vcf=temp("results/filtered/all.{vartype}_recalibrated.vcf"),
        idx=temp("results/filtered/all.{vartype}_recalibrated.vcf.idx"),
        tranches=temp("results/filtered/all.{vartype}_recalibrated.tranches"),
    params:
        mode=get_vartype_mode,
        resources={
            "ensembl": {"known": True, "training": False, "truth": False, "prior": 2.0}
        },
        annotation=[
            "MQ",
            "QD",
            "DP",
            "SB",
            "FS",
            "SOR",
            "ReadPosRankSum",
            "MQRankSum",
            "InbreedingCoeff",
        ],
        extra=config["params"]["gatk"]["VariantRecalibrator"],
    resources:
        mem_mb=16000,
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log",
    wrapper:
        "v5.1.0/bio/gatk/variantrecalibrator"


rule apply_recalibration:
    input:
        vcf="results/filtered/all.{vartype}.vcf.gz",
        ref="resources/chromosomes.fasta",
        tranches="results/filtered/all.{vartype}_recalibrated.tranches",
        recal="results/filtered/all.{vartype}_recalibrated.vcf",
    output:
        vcf="results/filtered/all.{vartype}.recalibrated.vcf.gz",
        idx="results/filtered/all.{vartype}.recalibrated.vcf.gz.tbi",
    params:
        mode=get_vartype_mode,
        extra=config["params"]["gatk"]["ApplyVQSR"],
    resources:
        mem_mb=16000,
    log:
        "logs/gatk/applyvqsr/{vartype}.log",
    wrapper:
        "v5.1.0/bio/gatk/applyvqsr"


rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype=(
                "recalibrated" if config["filtering"]["vqsr"] else "hardfiltered"
            ),
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    params:
        extra="--CREATE_INDEX",
    log:
        "logs/picard/merge-filtered.log",
    wrapper:
        "v5.1.0/bio/picard/mergevcfs"
