## USING POPULATION ALLELE FREQUENCIES FOR COMPUTING INDIVIDUAL ADMIXTURE ESTIMATES: 

Inference of ancestry is an important aspect of disease association studies as well as for understanding population history. We have developed a fast and accurate method for estimating the admixture proportions for an individual's ancestry using genotype or sequence data and population allele frequencies from a set of parental/reference populations. The method can work with genotype data or sequence data (aligned sequence reads in a BAM file) derived from low-coverage whole-genome sequencing, exome-sequencing or even targeted sequencing experiments. The method uses the L-BFGS-B code  (a limited memory BFGS algorithm with bound constraints) for optimizing the likelihood function and is extremely fast. 

The method is described in the paper: "Fast individual ancestry inference from DNA sequence data leveraging allele frequencies from multiple populations". Vikas Bansal and Ondrej Libiger. published in BMC Bioinformatics 2015. 

## INPUT: 

1. sorted BAM file for sequence data or simple genotype file (rsid genotype pairs)
2. population allele frequencies for common SNPs (generated using HapMap3 genotypes or other genotype datasets) 

## OUTPUT:  

admixture coefficients for each reference population 



## HOW TO RUN THE PROGRAM:

python runancestry.py gives all options for the program 

1. rename Makefile.simple to Makefile 
2. run 'make all' 


3. analyzing bam file for ancestry: python runancestry.py -f hapmap3.10clusters.admixture.AF.1-22only --bam sample.sorted.bam -o output.prefix -p 2 


## NOTES

1. To run on bam files, you will need the binary to calculate genotype likelihoods that is part of the iAdmix package. The source code for this binary is not part of the iAdmix code yet. 
2. It is not recommended to run the program directly from a VCF since VCFs typically don't have information about reference genotypes (0/0) and this may bias the ancestry inference. 
