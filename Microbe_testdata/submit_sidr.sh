#!/bin/bash
#SBATCH -p highmem
#SBATCH --qos jlfierst
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=UA44_SIDR
#SBATCH -e %J.err
#SBATCH -o %J.out

/jlf/acmccormack1/SIDR3/SIDR2.0 -f DSM44188.genome.nextpolish.fasta -t /jlf/acmccormack1/Pipelines/Metagenome_Decontamination/Data/NCBI_Taxdump/ -b DSM44188.blast.out -c actinobacteria  -n 1 -a DSM44188.bam
