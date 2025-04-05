# Genomic analysis of the Ixworth chicken

To reproduce and transparently follow the results of the paper, Snakemake pipelines, including R scripts, are provided here.

If you are interested in the Ixworth chicken, please also refer to [a previous study](https://doi.org/10.1080/00071668.2023.2246142).

In accordance with the paper, whole-genome resequencing (\~26x) of 49 Ixworth chicken liver samples have been deposited to ENA with accession number [PRJEB89160](https://www.ebi.ac.uk/ena/browser/view/PRJEB89160).

## General

The workflow is divided into three Snakemake pipelines:

1.  Data retrieval of ENA projects [PRJEB30270](https://www.ebi.ac.uk/ena/browser/view/PRJEB30270) and [PRJEB89160](https://www.ebi.ac.uk/ena/browser/view/PRJEB89160) ([fetch](./fetch/README.md))
2.  [GATK](https://gatk.broadinstitute.org/hc/en-us)-based joint variant calling ([variant_calling](./variant_calling/README.md))
3.  Population genetics based statistical evaluation of the data ([popgen](./popgen/README.md))

Program versions are pinned using [conda](https://anaconda.org/anaconda/conda).

To follow the described steps, a (Linux) compute server or a cluster is recommended with at least 64 available CPU cores, 256GB RAM and 10TB hard disk space.

Smaller configurations could possibly work at the expense of compute time.

## Additional resources

Most of the raw output files are provided inside the `output` directory. Some intermediary files were omitted if file sizes were too large.

Note that these files are unstructured and were not refined for publication.

## Installation

In order to run the different pipelines, python3, snakemake, conda (or any conda resolver like mamba), and apptainer (previously singularity) are required.

To get started, clone this repository:

``` bash
git clone https://github.com/Kirschluft/ixworth_popgen.git
cd ixworth_popgen
```

As recommended, install `conda` on your machine (see [official guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)) and follow the instructions.

Further, install `apptainer` (see [official guide](https://apptainer.org/docs/admin/main/installation.html#install-unprivileged-from-pre-built-binaries)) on your machine.

To run the Snakemake pipelines, [`snakemake`](https://pypi.org/project/snakemake/), [`snakemake-wrapper-utils`](https://pypi.org/project/snakemake-wrapper-utils/) and [`pandas`](https://pypi.org/project/pandas/) are required. The preferred installation is using `pip` inside a `conda` environment:

``` bash
conda create -n smk -y python=3.11
conda activate smk
pip install snakemake==9.1.1 snakemake-wrapper-utils==0.7.2 pandas==2.2.3
snakemake --version # check installation
```

Finally, two `apptainer` containers have to be built (required for `popgen`, not provided here):

``` bash
# building the containers will take a while
cur_dir=$(pwd)
cd popgen/workflow/envs
apptainer build rehh.sif Singularity_rehh.def
apptainer build ne.sif Singularity_ne.def
cd $cur_dir
```

Afterwards, you should be ready to run all pipelines (inside the `conda` environment `smk`).

## Reproduction

The number of cores used in the pipelines can be adjusted according to the needs. RAM restrictions are (typically) not added to the pipeline. Therefore, $> 256GB$ of RAM are recommended.

To check any of the pipelines for issues before running them, please use the following command in the respective directories:

``` bash
snakemake -np
```

For slightly more in-depth descriptions of the different pipelines, please refer to the corresponding `README.md` files ([fetch](./popgen/README.md), [variant_calling](./variant_calling/README.md), and [popgen](./popgen/README.md)).

### 1. Data retrieval

To retrieve the raw reads (`fastq`) from ENA, run the following command from the `fetch/` directory:

``` bash
snakemake --software-deployment-method conda apptainer --cores 1 --keep-going
```

This will download all ENA files into the `fetch/results/` directory using [`enabrowsertools`](https://github.com/enasequence/enaBrowserTools).

The total file sizes should be about 2.4TB and therefore, downloading will take a while. Note that some runs from project `PRJEB30270` are not downloaded since they were unavailable at the time of the study.

Once finished without errors, the files should be moved to the directory `variant_calling/data` for the variant calling pipeline. To achieve this, run the following command from the `root` directory (file extensions should be `_R1.fastq.gz` and `_R2.fastq.gz`):

``` bash
mkdir -p variant_calling/data/
for f in $(find fetch/results -type f -regex ".*_R?[12]\.\(fastq\|fq\)\.gz"); do
    base=$(basename "$f")
    newname=${base/_1.fastq.gz/_R1.fastq.gz}
    newname=${newname/_2.fastq.gz/_R2.fastq.gz}
    mv "$(realpath "$f")" "variant_calling/data/$newname"
done
```

### 2. Variant calling

To call variants using the GATK-based pipeline, run the following command from the `variant_calling/` directory:

``` bash
snakemake --software-deployment-method conda apptainer --cores 64 --scheduler greedy
```

This will generate an output directory at `variant_calling/results` that contains the final [`VEP`](https://mart.ensembl.org/info/docs/tools/vep/index.html) annotated (unfiltered) `VCF`-file (`variant_calling/results/final/all.ann.vcf.gz`) as well as a directory containing a quality control overview (`variant_calling/results/qc/multiqc.html`) and some general statistics on the variants (`variant_calling/results/stats/`).

Depending on the used configuration (local or cluster), the pipeline may take more than a week to finish due to the various steps (trimming, quality control, alignment, variant calling, merging). Once finished, the resulting `VCF`-file should be moved to `popgen/data/results.vcf.gz` for the statistical analysis pipeline:

``` bash
mv variant_calling/results/final/all.ann.vcf.gz popgen/data/results.vcf.gz
```

### 3. Population genetics

The final step of the workflow includes a pipeline to conduct various statistical analyses (genetic diversity, differentiation, signals of selection, etc.) regarding population genetics using `VCFtools`, `BCFtools`, `plink`, `beagle` and various R-scripts and R-packages.

Download a `bed.gz`-file containing QTLs of chicken for the reference genome `GRCg6a` from [Animal QTLdb](https://www.animalgenome.org/cgi-bin/QTLdb/GG/download?d=jiJHCuBaSwTzUn57F6we5&c=251004%205&dnm=chicken%20%28bed%20format%29) (release #55 ( Dec 23, 2024)) and save the file to `popgen/data/Animal_QTLdb_release55_chickenGRCg6a.bed.gz`.

To run the different steps of the population genetics analysis, use the following command from the `popgen/` directory:

``` bash
snakemake --software-deployment-method conda apptainer --cores 64
```

The final results will be put into `popgen/results` and contain various output files such as intermediate results, summary statistics, tables (e.g. selective sweeps), and images (distribution/manhattan plots).
