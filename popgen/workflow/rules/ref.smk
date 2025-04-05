rule get_genes:
    output:
        vcf="results/resources/genes.gff3.gz",
    log:
        "logs/ref/get_genes.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all",
    cache: True
    wrapper:
        "v5.1.0/bio/reference/ensembl-annotation"
