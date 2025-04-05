
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

The pipeline will jointly call all samples that are defined (and merge different lanes/files of the same sample), following GATK best practices.

## Sample and unit sheet

* Add samples to `config/samples.tsv`. Only the column `sample` is mandatory, but any additional columns can be added.
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. For each unit, define platform (typically one of `ILLUMINA`, `SOLID`, `LS454`, `HELICOS` and `PACBIO`), and either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system).

## Configuration file

The `config.yaml` file contains parameters that can be adjusted for the pipeline:

- `samples` and `units`: paths to `samples.tsv` and `units.tsv` files. Does not need to be changed.
- `ref`: contains species identifier, release and genome build from [Ensembl](https://www.ensembl.org). This is relevant for retrieving the reference genome, gene annotations and the VEP cache. Currently, no custom genomes are supported.
- `filtering`: contains GATK-filtering mechanisms. Variant recalibration is not yet fully supported (false). Typical GATK-hard-filtering for SNPs and indels are provided. If variant recalibration should be used, refer to [VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036351392-VariantRecalibrator) and the [snakemake wrapper](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/gatk/variantrecalibrator/test/Snakefile) and adjust `params` for the VariantRecalibrator as well as the file under `workflow/rules/filtering.smk` accordingly.
- `processing`: contains a parameter for removing or just marking duplicates (`remove-duplicates`) and a potential restriction to variants in a certain region provided by a bed file (using `restrict-regions`).
- `params`: contains general parameters of programs of the pipeline (GATK calls and fastp quality filtering/trimming). 
- `vep`: contains plugins and additional parameters for VEP (see [plugins](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) and [parameters](https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html)).
