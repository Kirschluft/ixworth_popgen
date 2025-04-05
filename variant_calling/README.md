# GATK Variant Calling

Adopted fork of [dna-seq-gatk-variant-calling](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling) (DOI: https://doi.org/10.5281/zenodo.4677629).

# Snakemake workflow

This Snakemake pipeline implements the [GATK best-practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) for calling small germline variants.

## Usage

Put your data into the `data` folder and fill in the files provided in the `config` directory (`config.yaml`, `samples.tsv`, `units.tsv`) in accordance to your needs. For more information on the configuration, check out [README.md](./config/README.md) inside `config/`.

Some convenience programs for creating `config/samples.tsv` and `config/units.tsv` from the files inside `data` are provided in `scripts`.

Once paired-end reads are put into the `data/` directory, `units.tv` and `samples.tsv` can be created from the `root` directory with the following commands:

```bash
chmod +x scripts/units.sh scripts/sample.sh
./scripts/units.sh
./scripts/sample.sh
```

The scripts assume that files have `_R1.fastq.gz` and `_R2.fastq.gz` file endings and are located inside `data/`, and sample names are extracted from the basenames of the files.

### Testing

To test whether your configuration files and samples are properly set up, run the following:

```bash
snakemake -np
```

### Running

Once the test run success, you can run the pipeline using the following command (using conda environments inside apptainer):

```bash
snakemake --software-deployment-method conda apptainer --cores 64
```

- `cores`: number of threads to use for all jobs in total
- `--software-deployment-method conda apptainer`: programs are installed into an apptainer container using conda.

Once the pipeline finishes, results are usually provided in `results` with the final annotated vcf file in `results/final/all.ann.vcf.gz`.
The MultiQC output `results/qc/multiqc.html` provides summary statistics on quality control and `results/stats` contains some statistics on the called variants based on bcftools and picard.

### Debugging (if something goes wrong)

The pipeline can still partially fail if any program runs out of memory or due to other issues.
Check out the logs inside the `logs` directory if any part of the pipeline fails.
Once you addressed the issue, rerun snakemake with `--rerun-incomplete`.


Note: During parts of the pipeline, files will be split based on available chromosomes.
Other scaffolds are omitted.
This might potentially create a lot of jobs for snakemake which might be slow to resolve. 
If this becomes a problem or jobs timeout, try changing snakemake parameters such as 
`--group-components`, `--batch`, `--scheduler`, `--keep-going` (see also [Snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html)).
