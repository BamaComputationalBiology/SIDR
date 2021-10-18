## README ##

SIDR3 is located in //jlf/acmccormack1/SIDR3/

within the /jlf/acmccormack1/SIDR3/src/ directory you can find the python files:
analysis_new.py
analysis_old.py
analysis.py 


SIDR3 calls the analysis.py file when running. It does not need to be recompiled to run the new copy of analysis.py. 
Whatever version of the analysis that you want to run should be put here. 




##analysis_new.py 
contains a working version of SIDR3
this version takes in a fasta, NCBI_taxdump, blast results, and the reads aligned to the genome (.bam). You also tell it the target and nodes

options:
-f contigs.fasta
-t NCBI_Taxdump/
-b blast.out
-c nematoda ##target
-n 3 ## nodes?
-a reads.bam ##bwa alignment of reads


this takes unclassified reads and classifies them using the desicion tree function. 
output of the decision tree is written to tokeep.txt and toremove.txt 

however it doesnt print out the target and nottarget contigs. 
I need to change it so it writes those files. and adds those to tokeep.txt and toremove.txt



### analysis_new_pea.py
NEWEST VERSION
working version of SIDR3 with usable output. can discuss. 

######### working script 
#!/bin/bash
#SBATCH -p highmem
#SBATCH --qos jlfierst
#SBATCH -n 3
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=UA44_SIDR
#SBATCH -e %J.err
#SBATCH -o %J.out

/jlf/acmccormack1/SIDR3/SIDR2.0 -f UA44.fasta -t /jlf/acmccormack1/Pipelines/Metagenome_Decontamination/Data/NCBI_Taxdump/ -b UA44.blast.out -c nematoda -n 3 -a UA44.bam

#########





##analysis_old.py  (copied with PEA comments on 6/14/2021)
does. not. work. 

current issue is with run_analysis() function on line 156. 

line 159 - comments. looks like Classifiers is just a list of options for classification?
line 161 - i changed line 162 from SVC to NuSVC. The previous version used a combination of paramaters from both functions that matched neither program. I fixed this. 
https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
https://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVC.html

line 170 - starting here the software uses "testdata". no idea where this came form or how it should be formatted. Looks like it was being tested somehow, and never got fixed. 
 testdata probably needs to be fed either the full list of contigs and blast hits, or the unclassified contigs? not sure. very unclear. 


Other changes:
line 4 changed to import NuSVC due to change at line 162. 
 I made some other changes in the import lines but where all small things that made it more functional. 

Looking at all of this now, are you supposed to tell SIDR which classifier to use when you call it?  









