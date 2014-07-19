USING POPULATION ALLELE FREQUENCIES FOR COMPUTING INDIVIDUAL ADMIXTURE ESTIMATES: 

Inference of ancestry is an important aspect of disease association studies as well as for understanding population history. We have developed a fast and accurate method for estimating the admixture proportions for an individual's ancestry using genotype or sequence data and population allele frequencies from a set of parental/reference populations. The method can work with genotype data or sequence data (aligned sequence reads in a BAM file) derived from low-coverage whole-genome sequencing, exome-sequencing or even targeted sequencing experiments. The method uses the L-BFGS-B code  (a limited memory BFGS algorithm with bound constraints) for optimizing the likelihood function and is extremely fast. The program (binaries compiled for Linux x86_64) can be downloaded from attachments below. 

The method is described in the paper: "Fast individual ancestry inference from DNA sequence data leveraging allele frequencies from multiple populations". Vikas Bansal and Ondrej Libiger. (submitted) 


INPUT: 

1. sorted BAM file for sequence data or simple genotype file (rsid genotype pairs)
2. population allele frequencies for common SNPs (generated using HapMap3 genotypes or other genotype datasets) 

OUTPUT:  admixture coefficients for each reference population 


More details coming soon....
