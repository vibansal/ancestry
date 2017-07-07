## iAdmix: USING POPULATION ALLELE FREQUENCIES FOR COMPUTING INDIVIDUAL ADMIXTURE ESTIMATES: 

Inference of ancestry is an important aspect of disease association studies as well as for understanding population history. We have developed a fast and accurate method for estimating the admixture proportions for an individual's ancestry using genotype or sequence data and population allele frequencies from a set of parental/reference populations. The method can work with genotype data or sequence data (aligned sequence reads in a BAM file) derived from low-coverage whole-genome sequencing, exome-sequencing or even targeted sequencing experiments. The method uses the L-BFGS-B code  (a limited memory BFGS algorithm with bound constraints) for optimizing the likelihood function and is extremely fast. 

The method is described in the paper: "Fast individual ancestry inference from DNA sequence data leveraging allele frequencies from multiple populations". Vikas Bansal and Ondrej Libiger. published in BMC Bioinformatics 2015. 

## INPUT: 

1. sorted BAM file for sequence data or simple genotype file (rsid genotype pairs)
2. population allele frequencies for common SNPs (generated using HapMap3 genotypes or other genotype datasets) 

## OUTPUT:  

admixture coefficients for each reference population 


## HOW TO RUN THE PROGRAM:

To compile the code: run 'make all' in the directory with the source code. This should create the executables 'ANCESTRY' and 'calculateGLL' (in the sub-directory parsebam). 

python runancestry.py gives all options for the program 


1.  analyzing bam file for ancestry: python runancestry.py -f populations.frequencies.txt --bam sample.sorted.bam -o sample.output --path path_directory_with_executable 

2. Example for genotype file: python runancestry.py --freq populations.frequencies.txt --geno sample.genotypes --out sample.ancestry 

3. Example for plink genotype file: python runancestry.py --freq populations.frequencies.txt --plink sample.genotypes --out sample.ancestry

For plink, the program will assume that the files sample.genotypes.ped and sample.genotypes.map exist



## NOTES

1. # the allele frequency file should be sorted by chromosome and position #

2. For running iAdmix, the path to the directory where the 'ANCESTRY' executable is located needs to be provided using the --path option to runancestry.py. This should be the directory where you downloaded the source code and compiled it. 

3. To run on bam files, you will need to calculate genotype likelihoods using the reads that overlap the variant sites. iAdmix provides a program called 'calculateGLL' for doing this. The  binary file (compiled on ubuntu x86\_64 platform) is available in the github repository. The source code has recently been added to the github repository and can be compiled with the 'make all' command. 
 
4. It is not recommended to run the program directly from a VCF since VCFs typically don't have information about reference genotypes (0/0) and this may bias the ancestry inference. 

5. Make sure that the chromosome names ('chr1' vs '1') and the reference genome version (hg18 vs hg19) in the BAM file match the allele frequency file. If your chromosome names have the 'chr' prefix, use the command line option "--addchr=True" for the runancestry.py script 

6. The '-c' and '-m' options are experimental and only for genotype data with multiple individuals 
