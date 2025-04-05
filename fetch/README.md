# Fetching ENA project 

This directory contains a brief Snakefile that allows downloading [ENA](https://www.ebi.ac.uk/ena/browser/home) projects using `enaBrowserTools`.

Meta data for projects from ENA is required to download individual runs using `enaDataGet`.

## Usage

To change the project(s) to download, simply adjust the `config.yaml` and add a dict of accessions (e.g. PRJEB30270) and the corresponding meta accession file from ENA (tsv).

To retrieve the data of the defined projects, run the following command (after changing `config.yaml`):

```bash
snakemake --software-deployment-method conda apptainer --cores 1 --keep-going
```

After the pipeline has finished, the `results` directory should contain directories for each project (accession number) that contain the runs.

Notes:
- Downloading all files will take a while (total file sizes are ~1.6TB)
- `cores 1` is set so only a single process downloads the data to avoid overloading IO and network bandwidth
- `-keep-going` makes sure that all files are tried to be downloaded even if some fail
- If one of the downloading steps fails, check out the respective log file and retry the command until it works (ENA downloads fail sometimes)

### Debugging (if something goes wrong)

The pipeline can still partially fail if any program runs out of memory or due to other issues.
Check out the logs inside the `logs` directory if any part of the pipeline fails.
Once you addressed the issue, rerun snakemake with `--rerun-incomplete`.
