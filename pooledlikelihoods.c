// first derivative of ancestry likelihood function...returns the partial derivatives at values x[] in array y[]
void ancestryCLLd_pooled(double x[],double deriv[])
{
	double ll=0.0, sum =0, WAFSUM=0.0; double WAFSUM_c = 0.0; double AFSUM = 0.0; double constLL =0;
	int i=0,j=0,k=0;
	for (j=0;j<POPS;j++) { AFSUM += x[j]; deriv[j] = -100000; }  
	constLL = (double)POOL_SIZE/AFSUM; constLL *= SNPS;
	double dsum = 0; double prior=0; double logvec[2];
	double logn[POOL_SIZE+1]; for (i=1;i<=POOL_SIZE;i++) logn[i] = log(i); double l1=0;

	for (i=0;i<SNPS;i++) 
	{
		WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j];  WAFSUM_c = AFSUM - WAFSUM;
		logvec[0] = log(WAFSUM); logvec[1] =log(WAFSUM_c);
		prior = logvec[0]*POOL_SIZE; ll = prior + GENOTYPES[i][0]; sum = ll; SCORES[0] = ll;
		for (j=1;j<=POOL_SIZE;j++) 
		{
			prior += logn[POOL_SIZE-j+1] - logn[j] - logvec[0] + logvec[1];
			ll = prior + GENOTYPES[i][j]; SCORES[j] =  ll; 
			if (sum > ll) sum += log(1.0+exp(ll-sum));  else sum = ll + log(1.0+exp(sum-ll));
		}

		for (k=0;k<POPS;k++)
		{
			ll = log(((double)POOL_SIZE*AFMATRIX[i][k])/WAFSUM) + SCORES[0]; dsum = ll; 
			for (j=1;j<=POOL_SIZE;j++)
			{ 
				l1 = (double)(POOL_SIZE-j)*AFMATRIX[i][k]*WAFSUM_c + (double)j*(1.0-AFMATRIX[i][k])*WAFSUM; 
				l1 /= WAFSUM*WAFSUM_c; 
				ll = log(l1) + SCORES[j]; 
				if (dsum > ll) dsum += log(1.0+exp(ll-dsum));  else dsum = ll + log(1.0+exp(dsum-ll));
			}
			dsum -= sum; 
			if (deriv[k] > dsum) deriv[k] += log(1.0+exp(dsum-deriv[k])); else deriv[k] = dsum + log(1.0+exp(deriv[k]-dsum));
		}
	}
	for (j=0;j<POPS;j++) 
	{
		deriv[j] = constLL - exp(deriv[j]); 
	}
	//for (j=0;j<POPS;j++) fprintf(stdout,"derivative %d %f %f \n",j,x[j],deriv[j]);
}

// instead of single value, take weighted sum of w_0.waf.waf + w_1.2.waf.(1-waf) + w_2.(1-waf).(1-waf) when we have genotype likelihoods
double ancestryCLL_pooled(double x[])
{
	double LL = 0.0,ll=0.0,sum=0; 
	int i=0,j=0; double WAFSUM=0.0; double WAFSUM_c = 0.0; double AFSUM = 0.0; double prior = 0; 
	for (j=0;j<POPS;j++) AFSUM += x[j];  
	double constLL = (double)POOL_SIZE*log(AFSUM); 
	double dsum = 0; double logvec[2];

	double logn[POOL_SIZE+1]; for (i=1;i<=POOL_SIZE;i++) logn[i] = log(i);

	for (i=0;i<SNPS;i++) 
	{
		WAFSUM =0; for (j=0;j<POPS;j++) WAFSUM += x[j]*AFMATRIX[i][j];  WAFSUM_c = AFSUM - WAFSUM;
		logvec[0] = log(WAFSUM); logvec[1] =log(WAFSUM_c);
		prior = logvec[0]*POOL_SIZE; ll = prior + GENOTYPES[i][0]; sum = ll; SCORES[0] = ll;
		for (j=1;j<=POOL_SIZE;j++) 
		{
			prior += logn[POOL_SIZE-j+1] - logn[j] - logvec[0] + logvec[1];
			ll = prior + GENOTYPES[i][j]; SCORES[j] =  ll; 
			if (sum > ll) sum += log(1.0+exp(ll-sum));  else sum = ll + log(1.0+exp(sum-ll));
		}
		LL += sum -constLL;
	}
	return -1*LL;
}
void ancestryCLLd1(double x[], double deriv[])
{
	double ll_orig = ancestryCLL_pooled(x); double ll_new = 0;
	int i=0,j=0; double delta = 1e-8; 
	for (j=0;j<POPS;j++)
	{
		x[j] += delta;
		ll_new = ancestryCLL_pooled(x); deriv[j] = ll_new-ll_orig; deriv[j] /= delta; 
		x[j] -= delta;	
	}	
}

