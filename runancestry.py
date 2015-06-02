#! /usr/bin/env python
# author Vikas Bansal vbansal@scripps.edu

## TODO 
## 1. print # of sites and average read-depth per site in calculateGLL code 
## 3. option to flip alleles for - strand SNPs
## 4. option to delete intermediate ancestry files created 
## 5. option to match SNPs using rsid or chr:position 
## 6. option to read VCF file with one or more individuals
## 7. option to read a file with list of BAM files for sequencing experiment 
## 8. option to specify --delta threshold for parsimonious ancestry estimation
## 9. along with admixture estimates: also print delta for each value...
## 10. read plink .bed files binary format for genotypes
## 11. option to do ancestry for only specific samples in ped file

## running code using plink files
## python ~/CODE/JOINTCODE-coral/ancestry/runancestry.py --plink simpar.AFR --freq ../HapMap-AF/hapmap3.255393.AFmatrix.b36 --out sim_AFR --path ~/CODE/JOINTCODE-coral/ancestry/

import sys, os, glob, string,math, time
from math import log
from optparse import OptionParser
from subprocess import call 
POOLSIZE= 2;
INCLUDE_STRAND_AMBIG = 0;
CORES =1;
WINDOW =0;
PARSIMONY =1;
HWE_CHECK = 0;
OUTPUT_GLL_ONLY = 0;

##################################################################################################################################

def output_GLL_inputfile(afmatrixfile,outfile):
	File1 = open(outfile,'w');
	File = open(afmatrixfile);
	for line in File: 
		if line[0] != '#':
			snp = line.strip().split(); 
			print >>File1, 'SNP',snp[0].strip('chr'),snp[1],snp[3],snp[4],snp[3] + '/' + snp[4],100; 
	File.close();
	File1.close();

def make_ancestry_inputfile_simple(afmatrixfile,GLLfile,outfile):
	File1 = open(outfile,'w');
	File = open(GLLfile);  GLLs = {};
	for line in File: snp= line.strip().split(); GLLs[(snp[POOLSIZE+1].strip('chr'),int(snp[POOLSIZE+2]))] = ' '.join(snp[0:POOLSIZE+1]);
	File.close();

	File = open(afmatrixfile);
	for line in File: 
		snp = line.strip().split(); 
		if line[0] == '#': 
			print >>File1, '#GLL',
			for i in xrange(5,len(snp)): print >>File1, snp[i],
			print >>File1, '\n',
		else:
			try: 
				ll =  GLLs[(snp[0].strip('chr'),int(snp[1]))]; 
				print >>File1, ll,
				for i in xrange(5,len(snp)): print >>File1, snp[i],
				print >>File1, '\n',
			except KeyError: pass; 
	File.close();
	File1.close();

## function that processed bam file for ancestry calculation
def ancestry_pipeline_bam(matrixfile,bamfile,outfile,executable_dir_path):

	## the python I/O can be removed in C 

	print >>sys.stderr, "making input file for getting genotype likelihoods"; 
	output_GLL_inputfile(matrixfile,outfile+'.forGLL');
	#python ../runancestry.py hapmap3.allchroms.shared.matrix.hg19  > hapmap3.allchroms.shared.matrix.hg19.forhapcut

	print >>sys.stderr, "calculating genotype likelihoods for bam file",bamfile,"poolsize is: ",POOLSIZE 
	if OUTPUT_GLL_ONLY ==1: 
		call([ executable_dir_path + "/calculateGLL --allsites 1 -p " + `POOLSIZE` + " --bam " + bamfile + " --variants " + outfile +".forGLL" + " --mmq 30 > " + outfile + ".GLL" ],shell=True)
	else: 
		call([ executable_dir_path + "/calculateGLL -p " + `POOLSIZE` + " --bam " + bamfile + " --variants " + outfile +".forGLL" + " --mmq 30 > " + outfile + ".GLL" ],shell=True)
		# ../../hapCUT/calculateGLL --bam DM00231.bwamem.MD.bam --variants hapmap3.allchroms.shared.matrix.hg19.forhapcut --mmq 30 > DM00231.bwamem.MD.GLL

		print >>sys.stderr, "making input file for ancestry calculations"; 
		make_ancestry_inputfile_simple(matrixfile,outfile + ".GLL",outfile + ".ancestry.input");
		#python runancestry.py ancestry-data/hapmap3.allchroms.shared.matrix.hg19 ancestry-data/DM00231.bwamem.MD.GLL > ancestry-data/DM00231.bwamem.MD.input
		print >>sys.stderr, "ancestry admixture calculations for",bamfile,"poolsize is ",POOLSIZE;
		call(["time " + executable_dir_path + "/ANCESTRY -p " + `POOLSIZE` + " --pr 1 -i " + outfile + ".ancestry.input > " + outfile + '.ancestry.out'],shell=True);
		#call(["rm -f " + bamfile + ".GLL"],shell=True);
		#call(["rm -f " + bamfile + ".forGLL"],shell=True);
		#call(["rm -f " + bamfile + ".ancestry.input"],shell=True);
		#./ANCESTRY ancestry-data/DM00231.bwamem.MD.input



##################################################################################################################################

## add support for plink ped files for multiple samples or processing a file with list of bamfile names
## function that reads in PLINK ped and map files and determines ancestry for each individual in ped file
def make_ancestry_inputfile_plink(pedfile,mapfile,afmatrixfile,outfile,executable_dir_path,modulo_CORE): 

	RC = {}; RC['A'] = 'T'; RC['T'] = 'A'; RC['C'] = 'G'; RC['G'] = 'C'; 
	gll_low = -4; gll_high = 0; gll_equal = -0.477121;

	## read the allele frequency matrix file into a hashtable once for all samples 
	print >>sys.stderr, "reading allele frequency file",afmatrixfile;
	File = open(afmatrixfile);lines =0; AFtable = {}; headerline = [];
	for line in File: 
		snp = line.strip().split(); 
		if line[0] == '#' or lines ==0: headerline = snp; 
		elif INCLUDE_STRAND_AMBIG ==1: AFtable[snp[2]] = snp; 
		elif snp[3] == 'A' and snp[4] == 'T': continue; 
		elif snp[3] == 'T' and snp[4] == 'A': continue; 
		elif snp[3] == 'C' and snp[4] == 'G': continue; 
		elif snp[3] == 'G' and snp[4] == 'C': continue; 
		else: AFtable[snp[2]] = snp; 
		lines +=1; 
	File.close()

	## read the map file once so that we have SNP rsid for each SNP
	print >>sys.stderr, "reading mapfile",mapfile;
	try: File = open(mapfile,'r'); 
	except IOError: print >>sys.stderr, "map file",mapfile,"not found"; return -1;
	SNPlist = [];
	for line in File:  
		snp = line.strip().split();  
		if snp[0][0] == '#': continue; ## BUG if first line has '# chrom', make sure the lines match up to ped file... 10/04/2014 
		SNPlist.append([snp[1],snp[0],snp[2]]);  
	File.close();
	####### 


	#### read the ped file line by line and calculate ancestry for each individual 
	File = open(pedfile,'r'); samples = 0;
	for line in File:
		if CORES > 1 and samples%CORES != modulo_CORE: samples +=1; continue; 
		GLLs = {};
		genotypes = line.strip().split(); 
		sampleid = genotypes[1]; outfilename = outfile+'.'+sampleid + '.input'; 
		File1 = open(outfilename,'w');

		# Family ID      Individual ID   FatherID MotherID Sex  phenotype #SIM_1 SIM_1 0 0 1 -9 C C C C A 

		print >>File1, '#GLL',
		for i in xrange(5,len(headerline)): print >>File1, headerline[i],
		print >>File1, '\n',

		overlapping_markers =0; missing_markers = 0;
		for i in xrange(6,len(genotypes),2):
			s = i-6; s /=2; 
			try: 
				snp = AFtable[SNPlist[s][0]]; flag = 0;
				if genotypes[i] == snp[3] and genotypes[i+1] == snp[3]: print >>File1, gll_high,gll_low,gll_low,
				elif genotypes[i] == snp[4] and genotypes[i+1] == snp[4]: print >>File1, gll_low,gll_low,gll_high,
				elif genotypes[i] == snp[3] and genotypes[i+1] == snp[4]: print >>File1, gll_low,gll_high,gll_low,
				elif genotypes[i+1] == snp[3] and genotypes[i] == snp[4]: print >>File1, gll_low,gll_high,gll_low,
				elif genotypes[i+1] == RC[snp[3]] and genotypes[i] == RC[snp[4]]: print >>File1, gll_low,gll_high,gll_low,
				elif genotypes[i] == RC[snp[3]] and genotypes[i+1] == RC[snp[4]]: print >>File1, gll_low,gll_high,gll_low,
				elif genotypes[i] == RC[snp[4]] and genotypes[i+1] == RC[snp[4]]: print >>File1, gll_low,gll_low,gll_high,
				elif genotypes[i] == RC[snp[3]] and genotypes[i+1] == RC[snp[3]]: print >>File1, gll_high,gll_low,gll_low,
				else: 
					flag = 1; missing_markers +=1;
					#print i,'missing',snp,SNPlist[s],genotypes[i],genotypes[i+1]
				if flag ==0: 
					for j in xrange(5,len(snp)): print >>File1, snp[j],
					print >>File1, '\n',
					overlapping_markers +=1;
			except KeyError: 
				#print >>sys.stderr, 'variant not found',SNPlist[s][0];
				pass; 
			
		File1.close();
		print >>sys.stderr, "\n\nancestry admixture calculations for individual:",samples+1,sampleid,'using',overlapping_markers,'markers';
		#print >>sys.stderr, "missing",missing_markers,len(genotypes)/2;

		if WINDOW ==0: call(["" + executable_dir_path + "/ANCESTRY -p 2 --pr " + `PARSIMONY` + " --HWE " + `HWE_CHECK` + " -i " + outfilename +  " > " + outfilename + ".ancestry"],shell=True);
		else: 
			print >>sys.stderr, "calling BFGS method in windows";
			call(["time " + executable_dir_path + "/ANCESTRY.1 -p 2 --pr " + `PARSIMONY` + " -i " + outfilename  + " > " + outfilename + ".ancestry"],shell=True);
		#call(["rm -f " + outfilename],shell=True);
		#call(["rm -f " + outfilename +  " " + outfilename + ".ancestry"],shell=True);
		samples +=1;
		#sys.exit();
		"""
		"""
	File.close();

####################################################################################################################################################


## genotype file is rsid AG pairs 	
def make_ancestry_inputfile_rsid(afmatrixfile,GLLfile,outfile):

	RC = {}; RC['A'] = 'T'; RC['T'] = 'A'; RC['C'] = 'G'; RC['G'] = 'C'; 
	#gll_low = -3.602; gll_high = -0.000021715; gll_equal = -0.477121255; ## assuming error rate = 1/20000 or 50 errors per million genotypes
	gll_low = -4; gll_high = -0.000086868; gll_equal = -0.477121;
	#gll_low = -9; gll_high = 0.00000; gll_equal = -0.477121;
	column =1;

	## read genotype file and store genotypes in hashtable 
	File1 = open(outfile,'w');
	GLLs = {};
	File = open(GLLfile);  
	for line in File: 
		snp= line.strip().split(); 
		if snp[column] == 'NN' or snp[column] == 'N/N' or snp[column] == 'NoCall' or snp[column] == '--': continue;
		if snp[column][1] == '/': genotype = snp[column][0] + snp[column][2]; 
		else: genotype = snp[column];
		GLLs[snp[0]] = genotype;
		
	File.close();

	lines = 0; valid_snps = 0; 
	File = open(afmatrixfile);
	for line in File: 
		snp = line.strip().split(); 
		if line[0] == '#' or lines ==0: 
			print >>File1, '#GLL',
			for i in xrange(5,len(snp)): print >>File1, snp[i],
			print >>File1, '\n',
		## ignore A/T and C/G snps that are strand ambiguous 
		elif INCLUDE_STRAND_AMBIG ==0 and snp[3] == 'A' and snp[4] == 'T': continue; 
		elif INCLUDE_STRAND_AMBIG ==0 and snp[3] == 'T' and snp[4] == 'A': continue; 
		elif INCLUDE_STRAND_AMBIG ==0 and snp[3] == 'C' and snp[4] == 'G': continue; 
		elif INCLUDE_STRAND_AMBIG ==0 and snp[3] == 'G' and snp[4] == 'C': continue; 
		elif snp[3] in RC and snp[4] in RC: ## the two alleles are [ACTG], ignore indel D/I alleles 
			try: 
				GENOTYPE = GLLs[snp[2]];  flag =0;
				if GENOTYPE[0] == snp[3] and GENOTYPE[1] == snp[3]: print >>File1, gll_high,gll_low,gll_low,
				elif GENOTYPE[0] == snp[4] and GENOTYPE[1] == snp[4]: print >>File1, gll_low,gll_low,gll_high,
				elif GENOTYPE[0] == snp[3] and GENOTYPE[1] == snp[4]: print >>File1, gll_low,gll_high,gll_low,
				elif GENOTYPE[1] == snp[3] and GENOTYPE[0] == snp[4]: print >>File1, gll_low,gll_high,gll_low,
				elif GENOTYPE[1] == RC[snp[3]] and GENOTYPE[0] == RC[snp[4]]: print >>File1, gll_low,gll_high,gll_low,
				elif GENOTYPE[0] == RC[snp[3]] and GENOTYPE[1] == RC[snp[4]]: print >>File1, gll_low,gll_high,gll_low,
				elif GENOTYPE[0] == RC[snp[4]] and GENOTYPE[1] == RC[snp[4]]: print >>File1, gll_low,gll_low,gll_high,
				elif GENOTYPE[0] == RC[snp[3]] and GENOTYPE[1] == RC[snp[3]]: print >>File1, gll_high,gll_low,gll_low,
				else: 
					#print >>File1, gll_equal,gll_equal,gll_equal,
					#print >>sys.stderr, GENOTYPE,snp;
					flag = 1;
				if flag ==0: 
					for i in xrange(5,len(snp)): print >>File1, snp[i],
					print >>File1, '\n',
					valid_snps +=1;
					#print >>sys.stderr, snp[2];
			except KeyError: pass; 
		lines +=1;
	print >>sys.stderr, "ancestry genotype file has",valid_snps,"markers";
	File.close();
	File1.close();

def ancestry_pipeline_genotypes(matrixfile,genotype_file,outfile,executable_dir_path):
	print >>sys.stderr, "making input file for ancestry calculations"; 
	make_ancestry_inputfile_rsid(matrixfile,genotype_file,genotype_file + ".ancestry.input");
	print >>sys.stderr, "ancestry admixture calculations",genotype_file
	if WINDOW ==0: call(["time " + executable_dir_path + "/ANCESTRY --LRT 7.68 -p 2 --pr 1 -i " + genotype_file + ".ancestry.input > " + outfile],shell=True);
	else: 
		print >>sys.stderr, "calling BFGS method in windows";
		call(["time " + executable_dir_path + "/ANCESTRY.1 -p 2 --pr 1 -i " + genotype_file + ".ancestry.input > " + outfile],shell=True);
	#call(["time " + executable_dir_path + "/ANCESTRY -p 2 --pr 1 -i " + genotype_file + ".ancestry.input > " + outfile],shell=True);
	#call(["time ./ANCESTRY --ea 1 -p 2 --pr 1 -i " + bamfile + ".ancestry.input > " + outfile],shell=True);
	#call(["rm -f " + bamfile + ".ancestry.input"],shell=True);


###########################################################################################################################



parser = OptionParser();
#parser.add_option("--type",dest="type",type="string",help="type of input file: bam/geno/genotype/VCF/ped",default="");
parser.add_option("-b","--bam",dest="BAMfile",type="string",help="input bam file name",default="");
parser.add_option("--geno",dest="genofile",type="string",help="input file name for simple genotype file: rsid AG on each line",default="");
parser.add_option("--vcf","--VCF",dest="vcffile",type="string",help="input file name for VCF file",default="");
parser.add_option("--plink",dest="pedfile",type="string",help="input file name for plink file in ped/map format (provide prefix of ped file name without the .ped)",default="");
parser.add_option("-f","--freq",dest="AFfile",type="string",help="allele frequency file",default="");
parser.add_option("-o","--out",dest="outfile",type="string",help="prefix of output file names",default="iADMIX.out");
parser.add_option("-p","--poolsize",dest="POOLSIZE",type="int",help="pool size for non-diploid samples, default = 2",default=2);
parser.add_option("--pr","--parsimony",dest="pr",type="int",help="parsimonious ancestry estimation 0/1, default = 1",default=1);
parser.add_option("--hwe","--HWE",dest="HWE",type="int",help="HWE deviation estimation 0/1, default = 0",default=0);
parser.add_option("--strand",dest="include_strand_ambiguous",type="int",help="<0/1> include SNPs that are strand ambiguous (A/T, C/G): default = 0 (not included)",default=0);
parser.add_option("--path",dest="path",type="string",help="path to directory with executables ANCESTRY and calculateGLL",default=".");
parser.add_option("-c","--cores",dest="CORES",type="int",help="number of cores for processing in parallel",default=1);
parser.add_option("-m","--modcore",dest="MOD_CORE",type="int",help="0/1/2/3..../CORES-1",default=0);
parser.add_option("-w","--windows",dest="WINDOW",type="int",help="0/1",default=0);
parser.add_option("-g","--gll",dest="GLL_ONLY",type="int",help="0/1",default=0); # only output genotypes likelihoods, no ancestry calculation...
(options,args) = parser.parse_args(); POOLSIZE = options.POOLSIZE;  CORES = options.CORES; HWE_CHECK = options.HWE; OUTPUT_GLL_ONLY = options.GLL_ONLY;
WINDOW = options.WINDOW;
INCLUDE_STRAND_AMBIG = options.include_strand_ambiguous
PARSIMONY = options.pr; 

#file = open(sys.argv[1],'rb'); while True: byte = file.read(1); print ord(byte); sys.exit();


if (options.BAMfile == "" and options.genofile == "" and options.pedfile == "") or options.AFfile == "": 

	print >>sys.stderr,"\n####program iAdmix for estimation of ancestry using population allele frequencies###\n";
	print >>sys.stderr, "python runancestry.py --bam sample.bam --geno sample.genotypes --freq allelefrequency_file --out outputfile --poolsize 2"; 
	print >>sys.stderr, "\nExample: python runancestry.py --freq hapmap3.allchroms.shared.matrix --geno HGDP01254.genotypes --out HGDP12054.ancestry";
	print >>sys.stderr, "Example: python runancestry.py --freq hapmap3.allchroms.shared.matrix --bam HGDP01254.sorted.bam --out HGDP12054.ancestry";
	print >>sys.stderr, "Example: python runancestry.py --freq hapmap3.allchroms.shared.matrix --plink HGDP01254.genotypes (ped/map) --out HGDP12054.ancestry\n\n";
	parser.print_help()
	sys.exit();
elif options.pedfile != "": 
	make_ancestry_inputfile_plink(options.pedfile+'.ped',options.pedfile+'.map',options.AFfile,options.outfile,options.path,options.MOD_CORE);
	
elif options.BAMfile != "": 
	if not os.path.exists(options.path + "/ANCESTRY"): print >>sys.stderr, "ANCESTRY executable does not exist, please specify correct path with the executable using --path option\n"; sys.exit(); 
	if not os.path.exists(options.path + "/calculateGLL"): print >>sys.stderr, "calculateGLL executable does not exist, please specify correct path with the executable using --path option\n"; sys.exit(); 
	ancestry_pipeline_bam(options.AFfile,options.BAMfile,options.outfile,options.path);
elif options.genofile != "": 
	if not os.path.exists(options.path + "/ANCESTRY"): print >>sys.stderr, "ANCESTRY executable does not exist, please specify correct path with the executable using --path option\n"; sys.exit(); 
	ancestry_pipeline_genotypes(options.AFfile,options.genofile,options.outfile,options.path);


########################################################################################################################################


