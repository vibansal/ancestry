/*
CODE for estimating ancestry of an individual sample (genotype likelihoods) using allele frequency for multiple populations 
author: Vikas Bansal, vbansal@scripps.edu
first implemented 12/01/12 
last modified 12/09/13
##chrom position rsid A1 A2 YRI CHB CHD TSI MKK MEX LWK CEU ASW JPT GIH

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Lbfgsb.3.0/lbfgsb.h"

//gcc Lbfgsb.3.0/lbfgsb.o -lm -o ANCESTRY ancestry.c
// we can speed up convergence by relaxing the stopping criteron

int REMOVE_POPS = 0; // if this variable is set to 0 only global likelihood maximization with all populations will be done
int POOL_SIZE = 2; 
int SNPS= 10; int POPS = 6; 
double** AFMATRIX; double** GENOTYPES; double* SCORES;
int SIMULATE = 0;
double LRT_threshold = 15.1367; // threshold for considering population significant in parsimonious ancestry fit
int INDIV_ID = 1;
int EXCLUDE_ADMIXED_POPS = 0;
int USE_DERIVATIVE = 0;
double BFGS_DELTA = 0.00001;
double MIN_ADMIXTURE = 0.001;
int RANDOM_START = 1;
int ADD_NOISE = 200;
int HWE_CHECK =0;
//double pvalues[3] = {0.01,0.001,0.0001}; double chivals[3] = {6.635,10.8275,15.1367};

#define MIN_AF 1e-3

#include "pooledlikelihoods.c" // separate likelihood functions for pooled genotypes

//store the admixture coefficients, genotypes, population allele frequencies in a single array x[]
// first is 'POPS' coefficients, 'SNPS' genotypes, followed by 'SNPS' frequencies for pop1, pop2 ......

// we could combine the two functions into a single one -> speedup

// first derivative of ancestry likelihood function...returns the partial derivatives at values x[] in array y[]
void ancestryCLLd(double x[],double deriv[])
{
	double ll=0.0, sum =0, WAFSUM=0.0; double WAFSUM_c = 0.0; double AFSUM = 0.0; double constLL =0;
	int i=0,j=0; 
	for (j=0;j<POPS;j++) { AFSUM += x[j]; deriv[j] = 0.0; } 
	constLL = 2.0/AFSUM; 

	for (i=0;i<SNPS;i++) 
	{
		WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j];  WAFSUM_c = AFSUM - WAFSUM; 
		sum = WAFSUM*WAFSUM*GENOTYPES[i][0] + 2*WAFSUM*WAFSUM_c*GENOTYPES[i][1] + WAFSUM_c*WAFSUM_c*GENOTYPES[i][2];
		if (HWE_CHECK ==1) sum += x[POPS]*WAFSUM*WAFSUM_c*(GENOTYPES[i][0] + GENOTYPES[i][2]-2.0*GENOTYPES[i][1]);

		for (j=0;j<POPS;j++)
		{
			ll = WAFSUM*AFMATRIX[i][j]*GENOTYPES[i][0];
			ll +=  (WAFSUM*(1.0-AFMATRIX[i][j]) + WAFSUM_c*AFMATRIX[i][j])*GENOTYPES[i][1];
			ll += WAFSUM_c*(1.0-AFMATRIX[i][j])*GENOTYPES[i][2]; 
			if (HWE_CHECK ==1) ll += 0.5*x[POPS]*(WAFSUM*(1.0-AFMATRIX[i][j]) + WAFSUM_c*AFMATRIX[i][j])*(GENOTYPES[i][0] + GENOTYPES[i][2]-2.0*GENOTYPES[i][1]);
			deriv[j] -= 2.0*ll/sum-constLL; 
			// constLL is the same for all SNPs and all pops so the derivative is not correct....
		}
		if (HWE_CHECK ==1)
		{
			if (i==0) deriv[POPS] = 0; 
			deriv[POPS] -= (WAFSUM*WAFSUM_c*(GENOTYPES[i][0] + GENOTYPES[i][2]-2.0*GENOTYPES[i][1]))/sum;
		}

	}
	//for (j=0;j<POPS;j++) fprintf(stdout,"derivative %d %f %f \n",j,x[j],deriv[j]);
}

// does not account for finite sample size to estimate allele frequencies, error in AF
// we don't do maximum likelihood -> sampling is worse than ML... 
// mean and variance of CLL for given admixture proportions...
void ancestryCLL_expected(double x[],double* CLL_mean,double* CLL_variance,int iters)
{
	double LL = 0.0,ll=0.0,sum=0; 
	int i=0,j=0,iter=0; double WAFSUM=0.0; double WAFSUM_c = 0.0; double AFSUM = 0.0; double mean =0,var=0;
	for (j=0;j<POPS;j++) AFSUM += x[j]; for (j=0;j<POPS;j++) x[j] /= AFSUM;
	*CLL_mean = 0; *CLL_variance = 0; 
	double w1=0,w2=0;


	// calculated empirical mean and variance and also shape of distribution...
	for (iter=0;iter<iters;iter++)
	{
		for (i=0;i<SNPS;i++) 
		{
			WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j]; WAFSUM_c = 1.0-WAFSUM; 
			//WAFSUM *= WAFSUM; WAFSUM_c *= WAFSUM_c;
			//sum = WAFSUM*WAFSUM + 4*WAFSUM*WAFSUM_c + WAFSUM_c*WAFSUM_c;
			mean = 2*WAFSUM*WAFSUM*log(WAFSUM) + 2*WAFSUM*WAFSUM_c*log(2*WAFSUM*WAFSUM_c) + 2*WAFSUM_c*WAFSUM_c*log(WAFSUM_c);

			var = WAFSUM*WAFSUM*(2*log(WAFSUM)-mean)*(2*log(WAFSUM)-mean);
			var += 2*WAFSUM*WAFSUM_c*(log(2*WAFSUM*WAFSUM_c)-mean)*(log(2*WAFSUM*WAFSUM_c)-mean);
			var += WAFSUM_c*WAFSUM_c*(2*log(WAFSUM_c)-mean)*(2*log(WAFSUM_c)-mean);
			*CLL_mean += mean;
			*CLL_variance += var; 
		}
	}

	int bins = 200; double** counts = calloc(sizeof(double*),bins); for (i=0;i<bins;i++) counts[i] = calloc(sizeof(double),6); 
	for (i=0;i<bins;i++) 
	{
		for (j=0;j<4;j++) counts[i][j] = 0.0;
	}

	// this code is only called if iters ==0
	for (i=0;i<SNPS && iters ==0;i++) 
	{
		WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j]; WAFSUM_c = 1.0-WAFSUM; 
		if (WAFSUM > 0.5) // collapse 0.3 and 0.7 into same bin
		{
			WAFSUM = 1.0-WAFSUM; 
			// swap genotype likelihoods as well
			ll = GENOTYPES[i][0]; GENOTYPES[i][0] = GENOTYPES[i][2]; GENOTYPES[i][2] = ll; 
		}

		counts[(int)(WAFSUM*bins)][0] +=1; 
		counts[(int)(WAFSUM*bins)][1] += 0.5*GENOTYPES[i][1] + GENOTYPES[i][0]; 
		counts[(int)(WAFSUM*bins)][2] += (0.5*GENOTYPES[i][1] + GENOTYPES[i][0])*(0.5*GENOTYPES[i][1] + GENOTYPES[i][0]);
		if (GENOTYPES[i][0] > 0.99) counts[(int)(WAFSUM*bins)][3] +=1;
		else  if (GENOTYPES[i][1] > 0.99) counts[(int)(WAFSUM*bins)][4] +=1;
		else if (GENOTYPES[i][2] > 0.99) counts[(int)(WAFSUM*bins)][5] +=1;

		//WAFSUM *= WAFSUM; WAFSUM_c *= WAFSUM_c;
		//sum = WAFSUM*WAFSUM + 4*WAFSUM*WAFSUM_c + WAFSUM_c*WAFSUM_c;
		mean = 2*WAFSUM*WAFSUM*log(WAFSUM) + 2*WAFSUM*WAFSUM_c*log(2*WAFSUM*WAFSUM_c) + 2*WAFSUM_c*WAFSUM_c*log(WAFSUM_c);
		var = WAFSUM*WAFSUM*(2*log(WAFSUM)-mean)*(2*log(WAFSUM)-mean);
		var += 2*WAFSUM*WAFSUM_c*(log(2*WAFSUM*WAFSUM_c)-mean)*(log(2*WAFSUM*WAFSUM_c)-mean);
		var += WAFSUM_c*WAFSUM_c*(2*log(WAFSUM_c)-mean)*(2*log(WAFSUM_c)-mean);
		*CLL_mean += mean;
		*CLL_variance += var; 
		//fprintf(stdout,"SNP_NO %d af %f LL %f sum %f %f\n",i,AFMATRIX[i][0],LL,log(sum),sum);
	}

	// how to calculate variance of Y (predicted variable) -> variance of each bin and take weighted average ??
	// why excess counts in 0.99 range ??
	// two PDFs (true and estimated..) for allele frequency bin counts...  discretized...
	// variance of Y (genotype) conditional on predictors (allele frequency linear combination) = residual... 

	double ssq =0,mean_af =0,var_af=0,af_bin,chi_sq =0;

	for (i=0;i<bins;i++)
	{
		if (counts[i][0] < 100) continue;
		af_bin = (double)i/bins; 
		mean_af = counts[i][1]/(counts[i][0]); var_af = counts[i][2]/(counts[i][0]) - mean_af*mean_af;
		ssq += counts[i][0]*(mean_af-af_bin)*(mean_af-af_bin); 
		chi_sq = 0;

		//fprintf(stderr,"bin %d %d %f %f %f ssq %f var_bin %f\n",i,(int)counts[i][0],counts[i][1],mean_af,af_bin,ssq,var_af);
                //fprintf(stderr,"bin %d freq %0.4f obs_freq %0.4f size %d counts: %d,%d,%d expected: %0.2f,%0.2f,%0.2f \n",i,af_bin,mean_af,(int)counts[i][0],(int)counts[i][3],(int)counts[i][4],(int)counts[i][5],counts[i][0]*af_bin*af_bin,counts[i][0]*2*af_bin*(1.0-af_bin),counts[i][0]*(1.0-af_bin)*(1.0-af_bin));
	}
	fprintf(stderr,"bins %d ssq %f\n",bins,ssq);

	for (i=0;i<bins;i++) free(counts[i]);  free(counts);
}

// instead of single value, take weighted sum of w_0.waf.waf + w_1.2.waf.(1-waf) + w_2.(1-waf).(1-waf) when we have genotype likelihoods
double ancestryCLL(double x[])
{
	double LL = 0.0,ll=0.0,sum=0; 
	int i=0,j=0; double WAFSUM=0.0; double WAFSUM_c = 0.0; double AFSUM = 0.0; double constLL =0;
	for (j=0;j<POPS;j++) AFSUM += x[j]; 
	//constLL = 2.0*log(AFSUM); // sum of admixture proportions, divide to ensure sum = 1
	for (i=0;i<SNPS;i++) 
	{
		WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j];  WAFSUM /= AFSUM; WAFSUM_c = 1.0-WAFSUM;
		//WAFSUM_c = AFSUM - WAFSUM;
		sum = WAFSUM*WAFSUM*GENOTYPES[i][0] + 2*WAFSUM*WAFSUM_c*GENOTYPES[i][1] + WAFSUM_c*WAFSUM_c*GENOTYPES[i][2];
		if (HWE_CHECK ==1) 
		{
			sum += x[POPS]*WAFSUM*WAFSUM_c*(GENOTYPES[i][0] + GENOTYPES[i][2]-2.0*GENOTYPES[i][1]);
			if (sum < 0) fprintf(stderr,"error sum is negative %f i %d F %f genotypes %f,%f,%f\n",sum,i,x[POPS],GENOTYPES[i][0],GENOTYPES[i][1],GENOTYPES[i][2]);
		}
		LL += log(sum);
	}
	//LL -= 2.0*SNPS*log(AFSUM);
	// we could add additional penalty function that penalizes small ancestry coefficients.... 
	// beta*Penalty where beta is ancestry coefficient...
	//for (j=0;j<POPS;j++){ //if (x[j] > 0.001) LL += 10*log(x[j]/AFSUM); }
	return -1*LL;
}

#define DELIM " \t\n"

void determine_matrix_size(char* afile,int* snps,int* populations)
{
	*snps =0; *populations = 0; 
	FILE* fp = fopen(afile,"r");
	char temp[10000];
	char delims[] = " "; char *result;
	while (fgets(temp,10000,fp) != NULL) 
	{
		if ((*snps) ==0)
		{
			// parse the line to determine number of populations in matrix
			//fprintf(stdout,"line %s ",temp);
			result = strtok(temp,DELIM); 
			while (result != NULL) 
			{
				//fprintf(stdout,"result is %s\n",result); 
				(*populations)++;
				result = strtok(NULL,DELIM);
			}
		}
		(*snps)++; 
	}
	fclose(fp);
	(*snps)--; (*populations)--;
	fprintf(stdout,"snps in file %d pops %d\n",*snps,*populations);	
}

// add noise to allele frequencies using binomial sampling of new frequency
void add_noise(double** allelefreq,int snps,int pops,int nsamples)
{
	int i=0,j=0,k=0; double newfreq = 0; 
	for (i=0;i<snps;i++)
        {
		for (j=0;j<pops;j++)
		{
			newfreq =0;
			for (k=0;k<nsamples;k++)
			{
				if (drand48() < allelefreq[i][j]) newfreq += 1;
			}
			allelefreq[i][j] = newfreq/nsamples; 
		}
	}
}

// if allele frequency is 0/1, use small estimate
void fix_zero_frequencies(double** allelefreq,int snps,int pops)
{
	int i=0,j=0,k=0; 
	for (i=0;i<snps;i++)
        {
		for (j=0;j<pops;j++) 
		{
			if (allelefreq[i][j] < MIN_AF) 
			{
				//fprintf(stderr,"err %f \n",allelefreq[i][j]);
				allelefreq[i][j] = MIN_AF; 
			}
			else if (1.0-allelefreq[i][j] < MIN_AF) allelefreq[i][j] = 1.0-MIN_AF; 
			//fprintf(stderr,"ddd \n");
		}
	}
}

void read_allelefreqs(char* afile,double** allelefreq,double** genotypes,char** poplabels)
{
	int i=0,j=0;
	FILE* fp = fopen(afile,"r");
	char temp[1024]; double log10_const = log(10.0);

	// scan first line for pop labels...
	fscanf(fp,"%s ",temp); for(j=0;j<POPS-1;j++) fscanf(fp,"%s ",poplabels[j]); fscanf(fp,"%s\n",poplabels[j]);

	for (i=0;i<SNPS;i++)
	{
		for(j=0;j<=POOL_SIZE;j++) fscanf(fp,"%lf ",&genotypes[i][j]);
		if (POOL_SIZE ==2)
		{
			if (SIMULATE ==0) 
			{ 
				for (j=0;j<3;j++) genotypes[i][j] = pow(10,genotypes[i][j]);
				for (j=0;j<3;j++) 
				{
					//if (genotypes[i][j] < 0.001) genotypes[i][j] = 0;  else if (genotypes[i][j] > 0.999) genotypes[i][j] = 1; 
				}
			} // assumed to be LL log base 10
		}
		else
		{
			for (j=0;j<=POOL_SIZE;j++) genotypes[i][j] *= log10_const;
		}
		for (j=0;j<POPS-1;j++) fscanf(fp,"%lf ",&allelefreq[i][j]);
		fscanf(fp,"%lf\n",&allelefreq[i][POPS-1]);
	}
	//if (ADD_NOISE > 0) add_noise(allelefreq,SNPS,POPS,ADD_NOISE);
	fix_zero_frequencies(allelefreq,SNPS,POPS);
	fclose(fp);
}

// in final solution, also output delta for removing this one population and replacing by others...
int find_parsimonious_admixture(double* fullvec,int* nbd,double* lowbound,double* upbound,char** poplabels,double* full_LL,double* pop_deltas)
{
        // chisquare of 10.8275 = 0.001 | 6.635 = 0.01 
        int iter=0; int numpars = POPS; int i=0; int non_zero_pops = 0;
        double temp_af=0; double sum=0; double maxval_global=0;
        double maxval[numpars]; double invec[numpars]; double solidpops[numpars]; int sp=0; double newcoeffs[numpars];
	double delta_LL=0;
	int negdelta = 0; // # of negative changes in likelihood in BFGS

        for (i=0;i<POPS;i++) { invec[i] = fullvec[i]; pop_deltas[i] = 0; } 
        // remove populations with low admixture coefficient    
        for (iter =0;iter < POPS;iter++)
        {
                if (invec[iter] < MIN_ADMIXTURE)
                {
                        nbd[iter] = 2; upbound[iter] = 0.0; invec[iter] = 0.0;
                }
                else non_zero_pops++;
                solidpops[iter] = 0;
        }
        double mindelta = 0; int pop_to_remove = -1; double adf=0; int it=0;
        while (non_zero_pops >= 2)
        {
                mindelta = 1000;
                for (iter =0;iter < POPS;iter++)
                {
                        if (upbound[iter] < MIN_ADMIXTURE || solidpops[iter] ==1) continue;
                        nbd[iter] = 2; upbound[iter] = 0.0; temp_af = fullvec[iter]; invec[iter] = 0.0;

			negdelta = 100; it =0;
			while (negdelta >=2 && ++it < 5) 
			{
				for (i=0;i<POPS;i++) 
				{
					if (upbound[iter] >= MIN_ADMIXTURE) invec[iter] = drand48();
				}
				if (POOL_SIZE ==2) maxval[iter] = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL,ancestryCLLd,lowbound,upbound,nbd,10,&negdelta);
				else if (USE_DERIVATIVE ==1) maxval[iter] = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,ancestryCLLd_pooled,lowbound,upbound,nbd,0,&negdelta);
				else maxval[iter] = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,NULL,lowbound,upbound,nbd,0,&negdelta);
			}

			delta_LL = *full_LL-maxval[iter];
			if (delta_LL < 0) 
			{
				fprintf(stderr,"negative delta %f suggests sub-optimal global optimization \n",delta_LL);
				//return -1;
			}
			pop_deltas[iter] = delta_LL;
                        fprintf(stdout,"maxval %f for POP %s delta %f ",maxval[iter],poplabels[iter],delta_LL);
                        sum=0; for (i=0;i<POPS;i++) sum += invec[i]; for (i=0;i<POPS;i++) invec[i] /=sum;
                        for (i=0;i<POPS;i++) fprintf(stdout,"%s:%0.4f ",poplabels[i],invec[i]); fprintf(stdout,"\n");
                        if (*full_LL-maxval[iter] < mindelta) 
			{ 
				mindelta = *full_LL-maxval[iter]; pop_to_remove = iter; adf = temp_af; 
                        	for (i=0;i<POPS;i++) newcoeffs[i] = invec[i];
			}
                        nbd[iter] = 1; upbound[iter] = 1.00; invec[iter] = temp_af;
                        if (*full_LL-maxval[iter] >= 10 && solidpops[iter]==0) { solidpops[iter] = 1;  sp++; }
                }
                if (sp ==0)
                {
                        fprintf(stderr,"difficult to estimate admixture coefficients with high confidence\n");
                        return 0;
                }

                if (mindelta*2 < LRT_threshold && upbound[pop_to_remove] > MIN_ADMIXTURE)
                {
                        non_zero_pops--;
                        fprintf(stderr,"contribution of ancestry from pop %s is not significant: delta %f admixture_prop %f\n",poplabels[pop_to_remove],mindelta,adf);
                        fprintf(stdout,"contribution of ancestry from pop %s is not significant: delta %f admixture_prop %f %d\n",poplabels[pop_to_remove],mindelta,adf,non_zero_pops);
                        iter = pop_to_remove; nbd[iter] = 2; lowbound[iter] = 0.0; upbound[iter] = 0.0; invec[iter] = 0.0;
                        *full_LL -= mindelta;
                        for (i=0;i<POPS;i++) fullvec[i] = newcoeffs[i];
                        for (i=0;i<POPS;i++) fprintf(stdout,"%s:%0.4f ",poplabels[i],newcoeffs[i]); fprintf(stdout,"\n");

                }
                else break;
        }
}


void simulate_admixed(double* coeffs,char* simparfile)
{
        int i=0,j=0; double sum=0,mean=0,sumsq=0,variance=0,maf;
        int allele0,allele1,genotype; char A0,A1;
	for (i=0;i<POPS;i++) coeffs[i] = 0; 
	
	FILE* fp = fopen(simparfile,"r");
	for(j=0;j<POPS-1;j++) fscanf(fp,"%lf ",&coeffs[j]); fscanf(fp,"%lf\n",&coeffs[j]); 
	fclose(fp);
	for (j=0;j<POPS;j++) fprintf(stdout,"c[%d]:%0.4f ",j,coeffs[j]); fprintf(stdout," admix-coeffs for simulation \n");
	//for (i=0;i<5;i++) coeffs[i] = 0.2;  for (i=0;i<10 && i < POPS-1;i++) coeffs[i] = 0.1;

	fprintf(stdout,"SIM_%d SIM_%d 0 0 1 -9",INDIV_ID,INDIV_ID);
        for (i=0;i<SNPS;i++)
        {
                maf=0;  for (j=0;j<POPS;j++) maf += coeffs[j]*AFMATRIX[i][j];
                if (drand48() < maf) allele0 = 0; else allele0 = 1;
                if (drand48() < maf) allele1 = 0; else allele1 = 1;
                genotype = allele0 + allele1;
		//fprintf(stdout,"genotypes %f %f ",GENOTYPES[i][0],GENOTYPES[i][1]);
		if ((int)GENOTYPES[i][0] == 0) A0 = 'A'; 
		else if ((int)GENOTYPES[i][0] == 1) A0 = 'C'; 
		else if ((int)GENOTYPES[i][0] == 2) A0 = 'G'; 
		else if ((int)GENOTYPES[i][0] == 3) A0 = 'T'; 
		if ((int)GENOTYPES[i][1] == 0) A1 = 'A'; 
		else if ((int)GENOTYPES[i][1] == 1) A1 = 'C'; 
		else if ((int)GENOTYPES[i][1] == 2) A1 = 'G'; 
		else if ((int)GENOTYPES[i][1] == 3) A1 = 'T'; 
		GENOTYPES[i][0] = 0;  GENOTYPES[i][1] = 0;  GENOTYPES[i][2] = 0; GENOTYPES[i][genotype] = 1; 
		if (allele0 == 0) fprintf(stdout," %c",A0);  else if (allele0 == 1) fprintf(stdout," %c",A1); 
		if (allele1 == 0) fprintf(stdout," %c",A0);  else if (allele1 == 1) fprintf(stdout," %c",A1); 
        }
	fprintf(stdout,"\n");

}


// if parsimony calculation indicates negative delta, restart from scratch

int main(int argc,char* argv[])
{
	time_t now; time(&now);    unsigned int iseed = (unsigned int)time(NULL);  srand48(iseed);
	int i=0,j=0; char* inputfile = NULL; char* simparfile = NULL;

        for (i=1;i<argc;i+=2)
        {
                if (strcmp(argv[i],"--input") ==0 || strcmp(argv[i],"-i") ==0) inputfile = argv[i+1];
                else if (strcmp(argv[i],"--snps") ==0)  SNPS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--POPS") ==0)  POPS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"-p") ==0)      POOL_SIZE = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--parsimony") ==0)      REMOVE_POPS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--pr") ==0)      REMOVE_POPS = atoi(argv[i+1]);
                else if (strcmp(argv[i],"--sim") ==0)    {   SIMULATE = 1; simparfile = argv[i+1]; } 
                else if (strcmp(argv[i],"--id") ==0)    {   INDIV_ID = atoi(argv[i+1]); } 
                else if (strcmp(argv[i],"--ea") ==0)    {   EXCLUDE_ADMIXED_POPS = atoi(argv[i+1]); } 
                else if (strcmp(argv[i],"--derivative") ==0)    {   USE_DERIVATIVE = atoi(argv[i+1]); } 
                else if (strcmp(argv[i],"--delta") ==0)    {   BFGS_DELTA = atof(argv[i+1]); } 
                else if (strcmp(argv[i],"--LRT") ==0)    {  LRT_threshold = atof(argv[i+1]); } 
                else if (strcmp(argv[i],"--random") ==0)    {  RANDOM_START = atoi(argv[i+1]); } 
                else if (strcmp(argv[i],"--HWE") ==0)    {  HWE_CHECK = atoi(argv[i+1]); } 
        }
        if (inputfile == NULL)
        {
                fprintf(stderr,"program requires at least one input file --input filename\n");
                return 1;
        }

        int s=0,p=0; determine_matrix_size(inputfile,&s,&p); SNPS = s; POPS = p;


	int numpars = POPS; if (HWE_CHECK ==1) numpars +=1;
	double lowbound[numpars]; double upbound[numpars];  int nbd[numpars]; 

	// we don't need to impose upper bound since we are anyway normalizing the coefficients | but this seems to speed up the algorithm a little bit
	for (i=0;i<POPS;i++) { nbd[i] = 1; lowbound[i] = 0.0; upbound[i] = 1.00; }

	double** allelefreq = calloc(SNPS,sizeof(double*)); for (i=0;i<SNPS;i++) allelefreq[i]=calloc(POPS,sizeof(double)); 
	double** genotypes =calloc(SNPS,sizeof(double*)); for (i=0;i<SNPS;i++) genotypes[i] = calloc(POOL_SIZE+1,sizeof(double));
	char** poplabels = calloc(POPS,sizeof(char*)); for (i=0;i<POPS;i++) poplabels[i] = calloc(1024,sizeof(char)); 

	read_allelefreqs(inputfile,allelefreq,genotypes,poplabels);

	int size_of_invec = POPS + SNPS + POPS*SNPS; 
	double *invec = calloc(POPS+1,sizeof(double)); 
	double *fullvec = calloc(POPS+1,sizeof(double)); double *pop_deltas = calloc(POPS+1,sizeof(double));

	double sum=0; double maxval_global=0; double F =0;
	int counter=0;
	AFMATRIX=allelefreq; GENOTYPES=genotypes;
	if (SIMULATE ==1) simulate_admixed(invec,simparfile); // simulate admixed individual 
	for (i=0;i<POPS;i++) invec[i] = (double)1.0/POPS; 

	if (HWE_CHECK ==1)
	{	
		lowbound[POPS] = 0.00; upbound[POPS] = 0.99; nbd[POPS] = 2; // inbreeding coefficient
		invec[POPS] = drand48();
	}

	//if (HWE_CHECK ==0) { lowbound[POPS]= 0; upbound[POPS] = 0; invec[POPS] = 0; } 
	//invec[0]=0.169; invec[1]=0.923; invec[2] = 0.259; invec[3] = 0.961; invec[4]=0.182; invec[5]=0.713; invec[6]=0.786; invec[7] = 0.679; invec[8]=0.751; invec[9]=0.137; invec[10]=0.938;
	SCORES = calloc(POOL_SIZE+1,sizeof(double));

	int negdelta =100; int it =0; double maxval_global_prev = 0;  int h=0;
	// rerun the method with new starting points if # of negative delta steps is 2 or more 
	// if value of likelihood in rerun is same as first iter, don't continue -> only maximum of two times typically...

	for (h=0;h<=HWE_CHECK;h++)
	{
		it =0; negdelta = 100; maxval_global_prev = 0; 
		if (h==1) { upbound[POPS] = 0.0; invec[POPS] = 0.00; numpars = POPS; }

		while (it < 5) 
		{
			for (i=0;i<POPS;i++) 
			{
				if (nbd[i] != 2 && (it > 0 || RANDOM_START ==1)) invec[i] = drand48(); // use random start if delta is negative
				else if (nbd[i] != 2 && it ==0) invec[i] = 1.0/POPS;
				//invec[i] =0; if (i==0 || i == POPS-1) invec[i] = 0.5; 
			}
			sum=0; for (i=0;i<POPS;i++) sum += invec[i]; for (i=0;i<POPS;i++) invec[i] /=sum;  
			maxval_global = ancestryCLL(invec);
			for (i=0;i<POPS;i++) fprintf(stdout,"%s:%0.3f ",poplabels[i],invec[i]); fprintf(stdout," initial sol %f\n",maxval_global);
			if (HWE_CHECK ==1 && h==0) fprintf(stdout,"initial Fval %f \n",invec[POPS]);

			if (POOL_SIZE ==2) maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL,ancestryCLLd,lowbound,upbound,nbd,10,&negdelta);
			else if (USE_DERIVATIVE ==1) maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,ancestryCLLd_pooled,lowbound,upbound,nbd,0,&negdelta);
			else maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,NULL,lowbound,upbound,nbd,10,&negdelta);
			it++;
			if (fabs(maxval_global-maxval_global_prev) < 0.001 || negdelta < 2) break;
			maxval_global_prev = maxval_global;
			//if (negdelta < 3) break;
		}
		if (HWE_CHECK ==1 && h==0) 
		{
			fprintf(stdout,"initial maxval %f ADMIX_PROP ",maxval_global);
			sum=0; for (i=0;i<POPS;i++) sum += invec[i]; for (i=0;i<POPS;i++) invec[i] /=sum;  
			for (i=0;i<POPS;i++) fprintf(stdout,"%s:%0.4f ",poplabels[i],invec[i]); fprintf(stdout,"\n");
			for (i=0;i<POPS;i++) { if (invec[i] >= MIN_ADMIXTURE) fprintf(stdout,"%s:%0.3f ",poplabels[i],invec[i]); } fprintf(stdout,"NON_ZERO\n");
			fprintf(stderr,"initial maxval %f ADMIX_PROP ",maxval_global); for (i=0;i<POPS;i++) fprintf(stderr,"%s:%0.4f ",poplabels[i],invec[i]); fprintf(stdout,"\n");
			fprintf(stderr,"\nFval %0.4f \n",invec[POPS]);
		}
	}

	//if (POOL_SIZE ==2) maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL,ancestryCLLd,lowbound,upbound,nbd,10); 
	//else maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,ancestryCLLd_pooled,lowbound,upbound,nbd,10);

	if (REMOVE_POPS ==1)
	{
		HWE_CHECK =0; numpars = POPS;
		find_parsimonious_admixture(invec,nbd,lowbound,upbound,poplabels,&maxval_global,pop_deltas);
		// final calculation of ancestry coefficients
		//if (POOL_SIZE ==2) maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL,ancestryCLLd,lowbound,upbound,nbd,10); 
		//else maxval_global = optim_LBFGS_BOUNDED(numpars,invec,ancestryCLL_pooled,ancestryCLLd_pooled,lowbound,upbound,nbd,10);
	}
	sum=0; for (i=0;i<POPS;i++) sum += invec[i]; for (i=0;i<POPS;i++) invec[i] /=sum;
        fprintf(stdout,"\nfinal maxval %f ADMIX_PROP ",maxval_global);
        for (i=0;i<POPS;i++) fprintf(stdout,"%s:%0.4f ",poplabels[i],invec[i]); fprintf(stdout,"\n");
        for (i=0;i<POPS;i++) if (invec[i] >= MIN_ADMIXTURE) { fprintf(stdout,"%s:%0.4f ",poplabels[i],invec[i]); }fprintf(stdout,"FINAL_NZ_PROPS\n");

        fprintf(stderr,"\nfinal maxval %f ADMIX_PROP ",maxval_global);
        for (i=0;i<POPS;i++) fprintf(stderr,"%s:%0.6f ",poplabels[i],invec[i]); fprintf(stderr,"\n");
        for (i=0;i<POPS;i++)  fprintf(stderr,"%s:%0.4f:%.2f ",poplabels[i],invec[i],pop_deltas[i]);  fprintf(stderr,"FINAL_ALL_PROPS\n");
        for (i=0;i<POPS;i++) if (invec[i] >= MIN_ADMIXTURE) { fprintf(stderr,"%s:%0.4f:%.2f ",poplabels[i],invec[i],pop_deltas[i]); }fprintf(stderr,"FINAL_NZ_PROPS ");


	//double AFSUM = 0; for (j=0;j<POPS;j++) AFSUM += invec[j]; 
	//for (j=0;j<POPS;j++) invec[j]= 0; invec[1] =1; //invec[5] = 0.2; 
	if (HWE_CHECK ==0)
	{
		double CLL_mean=0,CLL_variance=0; double explained=0;
		ancestryCLL_expected(invec,&CLL_mean,&CLL_variance,0); 
		explained = maxval_global-CLL_mean; explained /= SNPS; 
		if (explained > 0) explained = 0; 
		fprintf(stderr,"LL %f LL_exp %f %f %0.3f\n",maxval_global,CLL_mean,sqrt(CLL_variance),pow(10,explained));
	}
	else
	{
		F = invec[POPS]; invec[POPS] = F-0.05; if (invec[POPS] < 0) invec[POPS] = 0; 
		while (invec[POPS] < F + 0.05) 
		{
			maxval_global = ancestryCLL(invec); 
        		fprintf(stderr,"likelihood as function of F %f %f\n",invec[POPS],maxval_global);
			invec[POPS] += 0.005;
		}
	}


	// free memory	
	for (i=0;i<SNPS;i++) free(allelefreq[i]); for (i=0;i<SNPS;i++) free(genotypes[i]); free(allelefreq); free(genotypes);
	for (i=0;i<POPS;i++) free(poplabels[i]); free(poplabels); 
	free(invec); free(fullvec); free(SCORES);
		
	return(0);
}

