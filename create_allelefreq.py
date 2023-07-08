#! /usr/bin/env python
# author Vikas Bansal vbansal@scripps.edu

import sys, os, glob, string,math, time
from subprocess import call 

## read text file with allele frequency for each population, output allele-freq matrix 
def shared_allelefreqs(afile):
	File = open(afile);
	print >>sys.stderr, 'reading population allele frequencies';
	listofpopulations = {}; poptable = {}; snptable = {};
	lines = 0;
	for line in File:
		popsnp = line.strip().split();
		if popsnp[1] not in snptable: snptable[popsnp[1]] = [popsnp[2],int(popsnp[3]),popsnp[4],popsnp[7]]; 
		if popsnp[0] not in listofpopulations: listofpopulations[popsnp[0]] = 1; 
		poptable[(popsnp[1],popsnp[0])] = popsnp[5]; 
		lines +=1;
	File.close();
	print >>sys.stderr, 'reading population allele frequencies';

	AFmatrix = []; snps =0;
	for snp in snptable.iterkeys(): 
		missing = 0; AF = [];
		for pop in listofpopulations.iterkeys(): 
			try: AF.append(poptable[(snp,pop)]); 
			except KeyError: missing +=1;
		if missing ==0: 
			snpinfo = snptable[snp];
			AFmatrix.append([snpinfo[0],snpinfo[1],snp,snpinfo[2],snpinfo[3],AF]);
			snps +=1;


	AFmatrix.sort();
	print '#chrom','position','rsid','A1','A2',
	for pop in listofpopulations.iterkeys(): print pop,
	print;

	for i in xrange(snps):
		print AFmatrix[i][0],AFmatrix[i][1],AFmatrix[i][2],AFmatrix[i][3],AFmatrix[i][4],
		for f in AFmatrix[i][5]: print f,
		print;
	
	#print >>sys.stderr, 'noofpops',len(listofpopulations),'noofSNPs',lines/len(listofpopulations);
	return listofpopulations; # list of popID of all populations in the allele frequency file 


shared_allelefreqs(sys.argv[1]);
