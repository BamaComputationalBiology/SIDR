# SIDR
## README ##

SIDR3 is located on UAHPC in /jlf/acmccormack1/SIDR3/

This version takes in a fasta, NCBI_taxdump, blast results, and the reads aligned to the genome (.bam). You also tell it the target and nodes.
This takes unclassified reads and classifies them using the desicion tree function. 
output of the decision tree is written to tokeep.txt and toremove.txt 

Also prints fulldata.txt which includes the Sequence_Name, Label, GC_Content, and Coverage for each contig/sequence. 



### Example command: 
> SIDR2.0 -f contigs.fasta -t /path/to/NCBI_Taxdump/ -b blast.out -c nematoda -n 3 -a reads.bam 

### options:
> -f contigs.fasta ## your assembled genome
> 
> -t NCBI_Taxdump/  ## path to NCBI_Taxdump
> 
> -b blast.out  ## blast output file
> 
> -c nematoda ##target phylum (can optionally change the rank)
> 
> -n 3 ## nodes/threads??
> 
> -a reads.bam ## alignment of reads to assembled genome for coverage estimations

### optional options:
> -r rank (change default taxon rank from phylum)
> 
> -k kmers (change k-mer distribution from default)





### Python files
within the /src/ directory you can find the python files:
- analysis_new.py ## a working version not called by SIDR
- analysis_old.py ## an old version with other options, not working
- analysis.py ## the python script called by SIDR


SIDR3 calls the analysis.py file when running. It does not need to be recompiled to run the new copy of analysis.py. 
Whatever version of the analysis that you want to run should be put here. 


### working script for UAHPC from /jlf/peadams/SIDR3_runs/UA44/
> #!/bin/bash
> #SBATCH -p highmem
> #SBATCH --qos jlfierst
> #SBATCH -n 3
> #SBATCH --mem-per-cpu=8G
> #SBATCH --job-name=UA44_SIDR
> #SBATCH -e %J.err
> #SBATCH -o %J.out

> /jlf/acmccormack1/SIDR3/SIDR2.0 -f UA44.fasta -t /jlf/acmccormack1/Pipelines/Metagenome_Decontamination/Data/NCBI_Taxdump/ -b UA44.blast.out -c nematoda -n 3 -a UA44.bam



