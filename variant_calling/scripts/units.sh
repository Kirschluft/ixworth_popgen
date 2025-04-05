#!/bin/bash

printf "sample\tunit\tplatform\tfq1\tfq2\n" > ./config/units.tsv
declare -A unit_counts

for fq1 in data/*_R1.fastq.gz; do
    filename=$(basename "$fq1")
    sample=$(echo "$filename" | cut -d'_' -f1)
    underscores=${filename//[^_]/}
    count=${#underscores}

    if [[ "$count" -eq 1 ]]; then
        fq2="data/${sample}_R2.fastq.gz"
    else
        lane=$(echo "$filename" | cut -d'_' -f2)
        fq2="data/${sample}_${lane}_R2.fastq.gz"
    fi
    
    if [[ -f "$fq2" ]]; then
        unit_counts["$sample"]=$(( ${unit_counts["$sample"]} + 1 ))
        unit=${unit_counts["$sample"]}

        printf "${sample}\t${unit}\tILLUMINA\t${fq1}\t${fq2}\n" >> ./config/units.tsv
    else
        echo "Warning: Matching R2 file for $fq1 not found."
    fi
done
