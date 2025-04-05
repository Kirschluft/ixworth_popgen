
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

## Configuration file

The `config.yaml` file contains parameters that can be adjusted for the pipeline:

- `vcf_file`: path to the VCF-file that contains all populations/samples.
- `ref`: contains species identifier, release and genome build from [Ensembl](https://www.ensembl.org). This is relevant for retrieving the reference genome annotations. Currently, no custom genomes are supported.
- `populations`: contains a map of population identifiers and a corresponding `txt`-file that contains all sample names for that population (one per line), e.g. the name the first population is `ixworth` and its samples are in `data/populations/ixworth.txt`. Sample names should be the same as the ones used in the `VCF`-file.
- `n_chromosomes`: number of chromosomes (highest chromosome number).
- `chromosomes`: list of all chromosome or contig names that should be used within the pipeline. Chromosome/contig names should be the same as the ones used in the `VCF`-file. Using sex-chromosomes is currently not recommended since `beagle` haplotype estimation does not work with haploid (sex) regions.
- `windows`: window and step size (in bp) used by the various programs (for $F_{ST}$, $pi$, $D$, $iHS$, $XP-EHH$).
- `contrasts`: contrasts that should be used when estimating $F_{ST}$ (`fst`) or $XP-EHH$ (`haplo`). Should be lists of length two where both entries are population names defined in `populations`.
- `beagle`: additional parameters for running [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html). Currently, most parameters are left default, but some might be recommended in the future, e.g. recombination rate.
- `populations_map`: tsv file with two columns (Name, Population). Should map the sample identifiers of the `VCF`-file to population names.
- `genome_length`: total length of all the provided chromosomes. Can be calculated from faidx files.
- `prune`: plink filter parameters used for PCA/ADMIXTURE and contemporary $N_e$ calculations.
- `admixture`: parameters for the [Admixture](https://dalexander.github.io/admixture/) analysis.
- `qtls`: `bed.gz`-file containing QTLs from https://www.animalgenome.org/cgi-bin/QTLdb/index (here: chicken QTLs for GRCg6a)
- `ne`: additional parameters for historical $N_e$ estimations using smc++.
