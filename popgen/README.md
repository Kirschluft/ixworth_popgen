# Population genetics of Ixworth chicken

## General

This directory contains a Snakemake pipeline that runs various programs for data quality-filtering (VCF-preprocessing), structural analysis of populations (ADMIXTURE, PCA, decay of $LD$), summary statistics on populations, nucleotide diversity, inbreeding, effective Population size ($N_e$), and differentiation of populations ($pi$, $F_{ST}$) as well as indicators of selection ($F_{ST}$, $D$, $iHS$, $XP-EHH$).
The final results are indicators of (positive) selection and their related genes.

## Usage

Put your data into the `data` folder and fill in `config.yaml` provided in the `config` directory in accordance to your needs. For more information on the configuration, check out [README.md](./config/README.md) inside `config/`.

### Testing

To test whether your configuration files and samples are properly set up, run the following:

```bash
snakemake -np
```

### Running

Once the test run succeeds, you can run the pipeline using the following command:

```bash
snakemake --software-deployment-method conda apptainer --cores 64
```

- `cores`: number of threads to use for all jobs in total
- `--software-deployment-method conda apptainer`: programs are installed into an apptainer container using conda.

Once the pipeline finishes, results are provided in `results`.
Gathered summary statistics, images, etc. for each analysis are put in a corresponding folder (e.g. `pi`, `fst`, etc.) and typically named after their purpose.

### Debugging (if something goes wrong)

The pipeline can still partially fail if any program runs out of memory or due to other issues.
Check out the logs inside the `logs` directory if any part of the pipeline fails.
Once you addressed the issue, rerun snakemake with `--rerun-incomplete`.

## Scripts

Inside `workflow/scripts`, executable R-scripts (with a CLI) are provided for different parts of the analysis.
Theoretically, they can be used as stand-alone applications to conduct different statistical analyses (e.g. $F_{ST}$, $XP-EHH$, etc.), but usually require `workflow/util.R` (hard-coded path) to work.

## Further notes

- Ensembl is currently required to retrieve gene annotations.
- The input `VCF`-file is filtered according to GATK-based quality-filters, $< 50%$ missing genotypes, allele count of at least 2, the specified chromosomes/contigs, and only includes bi-allelic SNPs. The filtered data are used for all downstream analyses.
- QTL search can be integrated to compare regions (`bed`-file) of selection against known QTLs (e.g. [Animal QTLdb](https://www.animalgenome.org/cgi-bin/QTLdb/index)).
- For haplotype-based analyses, `Beagle` estimates haplotypes and also imputes the remaining missing genotypes. Currently, if alleles are missing, `bcftools +fixploidy` is used to fix them before running `Beagle`. Therefore, deep (high quality) WGS data is preferred.
- `iHS` analysis is conducted using the reference allele (here red junglefowl) as the ancestral allele.
- `XP-EHH` analysis is conducted using the first population as the reference (positive direction).
- LD-decay is estimated per chromosome and per population and smoothed using 1 kb bins.
- For ADMIXTURE and PCA analyses, SNPs are LD-pruned. ADMIXTURE might not support chromosome names that are non-numeric (like sex-chromosomes or contigs).
- Contemporary $N_e$ estimations are added with stricter linkage filters using NeEstimator (provided by dartRverse) on a few thousand markers.
- Historical $N_e$ estimations are performed using smc++ on all markers.
