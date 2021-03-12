/***********************************************************************
 *			     Ising Model 2D			       *
 *			     Pedro H Mendes 			       *
 *								       *
 *	 gcc -Wall ising_full.c -lgsl -lgslcblas -lm -static	       *
 *								       *
 **********************************************************************/

/***********************************************************************
 *                            INCLUDES                                 *
 **********************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>

/***********************************************************************
 *                            DEFINITIONS                              *
 **********************************************************************/
#define L	 16 
#define L2 	 (L*L)
#define TRAN	 100000 	//1e5
#define TMAX	 1000000	//1e6
#define J	 1.0

/***********************************************************************
 *                         GLOBAL VARIABLES                            *
 **********************************************************************/
int dE, M, ET;

/***********************************************************************
 *                            FUNCTIONS                                *
 **********************************************************************/
void initialize(gsl_rng *mt, double *boltz, int *spin, int *neigh, int *h_M, int *h_E, double TEMP);
void mc_routine(gsl_rng *mt, double *boltz, int *spin, int *neigh);
void print_histogram(FILE *file, int *hist, int SIZE, int STEP);

/***********************************************************************
 *                          MAIN PROGRAM                               *
 **********************************************************************/
int main(int argc, char *argv[])
{
	clock_t t_i, t_f;
	t_i = clock();

	int mcs;
	unsigned long int seed;
	char Arq1[100], Arq2[100], Arq3[100];
	FILE *arq1, *arq2, *arq3;
	
	double TEMP, CPU_TIME;
	TEMP = atof(argv[1]);
	
	int *spin, *neigh;
	int *h_M, *h_E;
	double *boltz;
	size_t size = L2*sizeof(int); 

	spin = (int*)malloc(size);
	neigh = (int*)malloc(size);
	h_M = (int*)malloc(2*size);
	h_E = (int*)malloc(4*size);
	boltz = (double*)malloc(sizeof(double)*9);

	seed = 0;

#ifdef RAND
	srand(time(NULL));
	seed = rand();
#endif
	gsl_rng *mt;
	mt = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(mt, seed);

	initialize(mt, boltz, spin, neigh, h_M, h_E, TEMP);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		mc_routine(mt, boltz, spin, neigh);
	}

	sprintf(Arq1, "temp_T%.2lfL%d.dsf", TEMP, L);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#seed = %ld\n#MCS\tM\tET\n", seed);

#ifdef DATA
	double s_e, s_e2, C;
	double s_m, s_m2, X;
	double med_m, med_e;

	s_e = 0.;
	s_e2 = 0.;
	C = 0.;
	s_m = 0.;
	s_m2 = 0.;
	X = 0.;
#endif

	for(mcs=0; mcs<TMAX; mcs++)
	{
		mc_routine(mt, boltz, spin, neigh);
	
		fprintf(arq1, "%d\t%d\t%d\n", mcs, M, ET);

		h_M[M + L2] += 1;
		h_E[ET + 2*L2] += 1;
#ifdef DATA
		s_e += (1.*ET);
		s_m += fabs(1.*M);
		s_e2 += (1.*ET)*(1.*ET);
		s_m2 += (1.*M)*(1.*M);
#endif
	}

#ifdef DATA
	s_e /= TMAX;
	s_m /= TMAX;

	med_m = s_m;
	med_e = s_e;

	s_e *= s_e;
	s_m *= s_m;
	s_e2 /= TMAX;
	s_m2 /= TMAX;
	
	C = (s_e2 - s_e)/(TEMP*TEMP);
	X = (s_m2 - s_m)/(TEMP);
#endif

	sprintf(Arq2, "hM_T%.2lfL%d.dsf", TEMP, L);
	sprintf(Arq3, "hE_T%.2lfL%d.dsf", TEMP, L);
	arq2 = fopen(Arq2, "w");
	arq3 = fopen(Arq3, "w");
	fprintf(arq2, "#seed = %ld\n#i\th_M[i]\n", seed);
	fprintf(arq3, "#seed = %ld\n#i\th_E[i]\n", seed);
	print_histogram(arq2, h_M, 2*L2, 2);
	print_histogram(arq3, h_E, 4*L2, 4);

	fclose(arq1);
	fclose(arq2);
	fclose(arq3);

	free(spin);
	free(neigh);
	free(h_M);
	free(h_E);
	free(boltz);

	gsl_rng_free(mt);

	t_f = clock();
	CPU_TIME = (double)(t_f - t_i)/CLOCKS_PER_SEC;

#ifdef DATA
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", TEMP, med_m, med_e, C, X, CPU_TIME);
#else
	printf("%lf\n", CPU_TIME);
#endif

	return 0;
}

/***********************************************************************
 *                          INITIALIZATION                             *
 **********************************************************************/
void initialize(gsl_rng *mt, double *boltz, int *spin, int *neigh, int *h_M, int *h_E, double TEMP)
{
	int i;

	boltz[4] = exp(-4.0*J/TEMP);
	boltz[8] = exp(-8.0*J/TEMP);

	for(i=0; i<L2; i++)
	{
		spin[i] = 2*(gsl_rng_uniform_int(mt,2)) - 1;
		neigh[i] = 0;
	}

	ET = 0.0;
	M = 0.0;

	for(i=0; i<L2; i++)
	{
		neigh[i] += spin[(i-L+L2)%L2];
		neigh[i] += spin[(i+1)%L + (i/L)*L];
		neigh[i] += spin[(i+L)%L2];
		neigh[i] += spin[(i-1+L)%L + (i/L)*L];

		ET += spin[i]*neigh[i];
		M += spin[i];
	}

	ET = (ET*(-J))/2.0;

	for(i=0; i<2*L2; i++)
	{
		h_M[i] = 0;
	}

	for(i=0; i<4*L2; i++)
	{
		h_E[i] = 0;
	}

	return;
}

/***********************************************************************
 *                   	   MONTE CARLO ROUTINE                         *
 **********************************************************************/
void mc_routine(gsl_rng *mt, double *boltz, int *spin, int *neigh)
{
	int i, t;
	double r, prob;
	
	for(t=0; t<L2; t++)
	{
		i = gsl_rng_uniform_int(mt, L2);
	
		dE = 0;

		dE = 2*spin[i]*neigh[i];

		if(dE <= 0)
		{
			spin[i] *= -1;
			ET = ET + dE;
			M = M + 2*spin[i];

			neigh[(i-L+L2)%L2] += 2*spin[i];
			neigh[(i+1)%L + (i/L)*L] += 2*spin[i];
			neigh[(i+L)%L2] += 2*spin[i];
			neigh[(i-1+L)%L + (i/L)*L] += 2*spin[i];
		}
		else
		{
			prob = boltz[dE];
			r = gsl_rng_uniform(mt);
		
			if(r < prob)
			{
				spin[i] *= -1;
				ET = ET + dE;
				M = M + 2*spin[i];
	
				neigh[(i-L+L2)%L2] += 2*spin[i];
				neigh[(i+1)%L + (i/L)*L] += 2*spin[i];
				neigh[(i+L)%L2] += 2*spin[i];
				neigh[(i-1+L)%L + (i/L)*L] += 2*spin[i];
			}
		}
	
	}

	return;
}

/***********************************************************************
 *                   	   PRINT HISTOGRAMS                            *
 **********************************************************************/
void print_histogram(FILE *file, int *hist, int SIZE, int STEP)
{
	int i;

	for(i=0; i<=SIZE; i+=STEP)
	{
		fprintf(file, "%d\t%d\n", (i-(SIZE/2)), hist[i]);
	}

	return;
}
