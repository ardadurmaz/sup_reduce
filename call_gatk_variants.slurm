#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --time=16:00:00
#SBATCH --mem=64gb

module load bwa
module load samtools
module load gatk/4.1.7.0
gatk HaplotypeCaller -I /home/axd497/wgs_ecoli/data/aligned/markdup_ecoli_merged.bam -O /home/axd497/wgs_ecoli/data/gatk_ecoli_calls.vcf -R /home/axd497/resources/ecoli_draft.fasta --sample-ploidy 1 --verbosity DEBUG
