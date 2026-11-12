#!/bin/bash
# run from the project root, after the sample metadata section of process_pig.R
set -e
sample_file="./custom_metadata/zhao_2021_sample_mapping.tsv"
fastq_dir="./raw_data/pig_atac_zhao/fastq"
aligned_dir="./processed_data/pig_atac_zhao/aligned"
log_dir="./processed_data/pig_atac_zhao/logs"
ref_index="./raw_data/susScr11/susScr11"
threads=8
mkdir -p "$fastq_dir" "$aligned_dir" "$log_dir"

#### Zhao et al. 2021 pig ATAC-seq: download, align to susScr11, delete fastq, one sample at a time ####
for srr in $(tail -n +2 "$sample_file" | cut -f1); do
        url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/0${srr: -2}/${srr}/${srr}"
        r1_fastq="${fastq_dir}/${srr}_1.fastq.gz"
        r2_fastq="${fastq_dir}/${srr}_2.fastq.gz"
        bam_file="${aligned_dir}/${srr}.sorted.bam"
        wget -c "${url}_1.fastq.gz" -O "$r1_fastq"
        wget -c "${url}_2.fastq.gz" -O "$r2_fastq"
        bowtie2 -p ${threads} -x "$ref_index" -1 "$r1_fastq" -2 "$r2_fastq" \
                --very-sensitive --no-discordant --no-mixed -X 2000 \
                2> "${log_dir}/${srr}.bowtie2.log" |
                samtools view -bS -q 30 -f 2 - |
                samtools sort -@ ${threads} -o "$bam_file" -
        samtools index "$bam_file"
        rm -f "$r1_fastq" "$r2_fastq"
done

#### Halstead et al. 2020 pig ATAC-seq peaks (GSE158414): union peaks as 200 bp windows for susScr11 to hg38 mapping ####
mkdir -p ./processed_data/pig_atac_halstead
zcat ./raw_data/pig_atac_halstead/peaks/*.bed.gz | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - |
        awk '{center = int(($2 + $3) / 2); start = center - 100; if (start < 0) start = 0; print $1"\t"start"\t"center + 100"\tpig_dhs_"NR}' \
        > ./processed_data/pig_atac_halstead/pig_dhs_200bp.bed
