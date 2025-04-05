import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")


report: "../report/workflow.rst"


container: "docker://continuumio/miniconda3:25.1.1-2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Helper functions #####


# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        df = pd.read_table(fai, header=None, usecols=[0], dtype=str, names=["ID"])
        return df.squeeze("columns")


# potentially causes issues if not two fastq files are given
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    return {"sample": [fastqs.fq1, fastqs.fq2]}


# Consider using units instead of sample for ID when different lane biases should be considered, etc.
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{unit}\tSM:{sample}\tPL:{platform}'".format(
        unit=wildcards.unit,
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    # paired-end sample
    return expand(
        "results/trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards
    )


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}.bam",
        sample=wildcards.sample,
    )


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_recal_input(bai=False):
    f = "results/dedup/{sample}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # need an index because random access is required
            f += ".bai"
            return f
        else:
            # no index needed
            return []
    else:
        return f


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    )


def get_vartype_mode(wildcards):
    return "SNP" if wildcards.vartype == "snvs" else "INDEL"


def get_filter(wildcards):
    return config["filtering"]["hard"][wildcards.vartype]


def get_read1(wildcards):
    """Get fastq files read 1 of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if "fq1" in fastqs.index:
        return [fastqs["fq1"]]
    else:
        return []


def get_read2(wildcards):
    """Get fastq files read 2 of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if "fq2" in fastqs.index:
        return [fastqs["fq2"]]
    else:
        return []


def get_units(wildcards):
    return units.loc[wildcards.sample].index
