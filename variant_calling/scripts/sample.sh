#!/bin/bash

echo "sample" > ./config/samples.tsv
for file in ./data/*.fastq.gz; do
    sample_name=$(echo $(basename "$file") | cut -d'_' -f1 | rev | cut -d'-' -f1 | rev)
    echo "$sample_name"
done | sort -u >> ./config/samples.tsv
