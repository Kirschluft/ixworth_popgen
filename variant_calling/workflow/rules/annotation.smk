rule annotate_variants:
    input:
        calls="results/filtered/all.vcf.gz",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
    output:
        calls=report(
            "results/final/all.ann.vcf.gz",
            caption="../report/vcf.rst",
            category="Calls",
        ),
        stats=report(
            "results/stats/all.stats.html",
            caption="../report/stats.rst",
            category="Calls",
        ),
    params:
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"],
    log:
        "logs/vep/annotate.log",
    threads: 4
    wrapper:
        "v5.9.0/bio/vep/annotate"  # v5.1.0 is no longer working
