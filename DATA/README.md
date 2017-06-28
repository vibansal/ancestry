

## The format of the allele frequency file required as input for iAdmix is as follows:


The first line should start with '#' and has information about the population names/ids. 

#chrom position rsid A1 A2 YRI CHB CHD TSI MKK LWK CEU JPT

Each subsequent line has information about the allele frequencies of the populations and the two alleles for each variant. Note that the allele frequency is the frequency of the 'A1' allele (column 4). It does not matter if it is the minor or major allele.


1 566875 rs2185539 T C 0.00343 0.00182 0.00229 0.00245 0.21240 0.00227 0.00223 0.00221

1 728951 rs11240767 T C 0.12410 0.00182 0.00231 0.00245 0.13820 0.11680 0.00223 0.00221

1 752721 rs3131972 G A 0.20550 0.76280 0.80280 0.84800 0.39100 0.34550 0.83480 0.72570



A sample allele frequency file is included in this folder. 
