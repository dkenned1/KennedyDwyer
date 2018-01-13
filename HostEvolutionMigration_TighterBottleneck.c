#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#define EPSILON 0.0001   // Define your own tolerance
#define DOUBLE_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define PI 3.141592653589793235028831971650288
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_sort.h>

struct Params
{
	double B;//=4.704513e-06;
	double phi;//=2.094968e-01;
	double c1;//=4.270473e+01;
	double c2;//=1.228317e+03;
	long int m;//=4.073841e+04;
	double sigmam;
	double sigmaN;
	double AppliedDose;
	double N;
};

/*
double DeathTime(int Genotype[], long int SNPs, double MeanKillTime, gsl_rng *r)
{
	return (gsl_ran_gamma(r, 3, 100));
}
*/
double DeathTime(double Genotypes[], long int Initial_Strains, double MeanKillTime, int Model, gsl_rng * r)
{
	struct Params Params;
	Params.B=4.704513e-06;
	Params.phi=2.094968e-01;
	Params.c1=4.270473e+01;
	Params.c2=1.228317e+03;
	Params.m=4.073841e+04;
//	Params.AppliedDose=844;
//	Params.AppliedDose=10000;
	Params.AppliedDose=156;
	Params.N=1e9;

	long long m=(long long)(Params.m);
	double N= Params.N;

	long int StrainOverflow=200;

	long int n;
	double MeanOfNextEvent;
	double NextEvent;
	double t=0;


	long int a;
	long int WhichVirus;
	long int Temp;

	long int Strains[StrainOverflow];

	double EffectiveDose=Params.c1*Params.AppliedDose/(Params.c2+Params.AppliedDose);

	t=0;
	n=gsl_ran_poisson (r, EffectiveDose);

	for (a=0; a<n; a++)
	{
		Strains[a]=1;
		if (Model==1)
		{
			Genotypes[a]=1.0/(double)n;
		}
	}
	for (a=n; a<StrainOverflow; a++)
	{
		Strains[a]=0;
		if (Model==1)
		{
			Genotypes[a]=0;
		}
	}

	while (1)
	{
		if (n==0)
		{
			return(1e100);
		}

		if (n==10000)
		{
			for (a=0; a<StrainOverflow; a++)
			{
				if (Model==2)
				{
					Genotypes[a]=(double)Strains[a]/10000.0;
				}
				else
				{

				}
			}
			return(t+log(N/(double)n)/Params.phi);
		}

		MeanOfNextEvent=1/(Params.phi*n+Params.B*n*m);
		t-=logf(gsl_rng_uniform_pos (r))*MeanOfNextEvent;
		NextEvent=gsl_rng_uniform_pos (r);

		WhichVirus=gsl_rng_uniform_int(r,n)+1;
		Temp=0;

		if (NextEvent<Params.phi*n*MeanOfNextEvent)
		{
			for (a=0; a<StrainOverflow; a++)
			{
				Temp+=Strains[a];
				if (Temp>=WhichVirus)
				{
					Strains[a]++;
					break;
				}
			}
			n++;
		}
		else
		{
			for (a=0; a<StrainOverflow; a++)
			{
				Temp+=Strains[a];
				if(Temp>=WhichVirus)
				{
					Strains[a]--;
					break;
				}
			}
			n--;
			m--;
		}
	}
}


long int Differences(int *Genotype1, int *Genotype2, long int SNPs)
{
	long int Polymorphisms=0;
	long int a;

	for(a=0; a<SNPs; a++)
	{
		if (Genotype1[a]!=Genotype2[a])
		{
			Polymorphisms++;
		}
	}

	return(Polymorphisms);
}

double fmin(double Value1, double Value2)
{
	if (Value1<Value2)
	{
		return(Value1);
	}
	else
	{
		return(Value2);
	}
}

void Recombination(int *P1, int *P2, int *F1, long int SNPs, gsl_rng *r)
{
	long int a;
	for (a=0; a<SNPs; a++)
	{
		if (gsl_rng_uniform(r)<.5)
		{
			F1[a]=P1[a];
		}
		else
		{
			F1[a]=P2[a];
		}
	}
}

void ShellSort(double a[], int l, int r)		//Sort large to small
{
	int i, j, k, h;
	double v;
    int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		13776, 4592, 1968, 861, 336,
		112, 48, 21, 7, 3, 1 };

    for ( k = 0; k < 16; k++)
    {
     	for (h = incs[k], i = l+h; i <= r; i++)
        {
          	v = a[i]; j = i;
          	while (j >= h && a[j-h] < v)
          	{
				a[j] = a[j-h]; j -= h;
			}
          	a[j] = v;
        }
	}
}


void OldShellSort(double a[], int l, int r)		//Sort small to large
{
	int i, j, k, h;
	double v;
    int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		13776, 4592, 1968, 861, 336,
		112, 48, 21, 7, 3, 1 };

    for ( k = 0; k < 16; k++)
    {
     	for (h = incs[k], i = l+h; i <= r; i++)
        {
          	v = a[i]; j = i;
          	while (j >= h && a[j-h] > v)
          	{
				a[j] = a[j-h]; j -= h;
			}
          	a[j] = v;
        }
	}
}

void WithinHostDrift(double TempHolder[], double InfectingStrain[], double StrainHolder[], long int Strains, int Model, gsl_rng *r)
{
	long int a, b;
	double ForceOfI, Temp;

	double StrainOverflow=200;

	if (Model==0)
	{
		for (a=0; a<Strains; a++)
		{
			TempHolder[a]=InfectingStrain[a];
		}

		for (b=0; b<StrainOverflow; b++)
		{
			Temp=gsl_rng_uniform(r);			//To sync random number generator up with other modesl
		}
	}
	else
	{
		for (a=0; a<Strains; a++)
		{
			TempHolder[a]=0;
		}

		for (b=0; b<StrainOverflow; b++)
		{
			Temp=gsl_rng_uniform(r);
			ForceOfI=0;
			for (a=0; a<Strains; a++)
			{
				ForceOfI+=InfectingStrain[a];
				if (ForceOfI>Temp)
				{
					TempHolder[a]+=StrainHolder[b];
					break;
				}
			}
		}
	}
}



int main(int argc, char *argv[])
{

	//printf("Here\n");
	long int StartTime=time(NULL);
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	long int Seed=atoi(argv[1]);
	gsl_rng_set(r, Seed);
	//printf("Seed=%ld\n", Seed);

	long int i,j,a,b,c;
	long int Years=100;
	double phi=atof(argv[3]);

	double gamma=0.2;
	double BaseReproduction=atof(argv[4]);

	double beta=.0002;					//5e-7;								//threshold is 1000 individuals
	double mu=0.017;	//Fuller et al. estimates mu=.4 per day, value comes from dividing this by 24 hours/day				//Old value was .1;
	double lambda=1e8;

	double ProbOfPredation;
	double TypeThree_a=0;			//0.967;
	double TypeThree_b=1000;			//0.14;

	double YearQuality;

	double CV=atof(argv[2]);
	long int NumberOfDeaths=1e6;

	long int Uninfected=-1;
	long int Dead=-10;

	int Model=atoi(argv[5]);

	long int Realizations=1;
	struct Params Params;
	Params.B=4.704513e-06;
	Params.phi=2.094968e-01;
	Params.c1=4.270473e+01;
	Params.c2=1.228317e+03;
	Params.m=4.073841e+04;
//	Params.AppliedDose=844;
	Params.AppliedDose=10000;
	Params.N=1e9;

	long int SNPs=50;
	long int Initial_Strains=SNPs;
	long int TempStrain;
	long int LargeNumber=1e7;

	long int TempIndex, TempOffspring;

	double **S=malloc(LargeNumber*sizeof(double *));	if (S==NULL){printf("Too large\n");}
	for (a=0; a<LargeNumber; a++)
	{
		S[a]=malloc(Initial_Strains*sizeof(double));	if (S[a]==NULL){printf("Too large\n");}
	}

	double **I=malloc(LargeNumber*sizeof(double*));	if (I==NULL){printf("Too large\n");}
	for (a=0; a<LargeNumber; a++)
	{
		I[a]=malloc(Initial_Strains*sizeof(double));	if (I[a]==NULL){printf("Too large\n");}
	}

	double **Old_I=malloc(LargeNumber*sizeof(double*));	if (Old_I==NULL){printf("Too large\n");}
	for (a=0; a<LargeNumber; a++)
	{
		Old_I[a]=malloc(Initial_Strains*sizeof(double));	if (Old_I[a]==NULL){printf("Too large\n");}
	}

	double **Total_V=malloc(LargeNumber*sizeof(double *)); if (Total_V==NULL){printf("Too large\n");}
	for (a=0; a<LargeNumber; a++)
	{
		Total_V[a]=malloc(Initial_Strains*sizeof(double)); if (Total_V[a]==NULL){printf("Too large\n");}
	}

	double MeanKillTime=1;

	double *DeathHolder=malloc(NumberOfDeaths*sizeof(double));	if (DeathHolder==NULL){printf("Too large\n");}

	double **StrainHolder=malloc(NumberOfDeaths*sizeof(double)); if (StrainHolder==NULL){printf("Too large\n");}
	for (a=0; a<NumberOfDeaths; a++)
	{
		StrainHolder[a]=malloc(200*sizeof(double)); if (StrainHolder[a]==NULL){printf("Too large\n");}
	}

	double *h=malloc(LargeNumber*sizeof(double));	if (h==NULL){printf("Too large\n");}

	unsigned int StrainNames[Initial_Strains];

	a=0;
	b=0;
	while (a<NumberOfDeaths)
	{
		DeathHolder[a]=DeathTime(StrainHolder[a], Initial_Strains, MeanKillTime, Model, r);

		if (DeathHolder[a]>10000)
		{
			if (b<(0.028*NumberOfDeaths))
			{
				b++;
				a++;
			}
			else
			{
			}
		}
		else
		{
			a++;
		}
		//printf("%lf\n", DeathHolder[a]);
	}


	double TempHolder[Initial_Strains];

	for (i=0; i<Realizations; i++)
	{

		//long int Strains=50;
		//long int MaxStrains=20000;

		long int TotalTotal_V=0;

		long int CadaverCounter=0;

		long int Initial_S=10000;
		long int Initial_I=5000;

		long int Total_E;


//		long int *Total_V=malloc(MaxStrains*sizeof(long int));	if (Total_V==NULL){printf("Too large\n");}
//		long int *Initial_Infectives=malloc(MaxStrains*sizeof(long int));	if (Initial_Infectives==NULL){printf("Too large\n");}

		//double *I=malloc(LargeNumber*sizeof(double));	if (I==NULL){printf("Too large\n");};



		long int Total_S=Initial_S;
		long int Total_I=Initial_I;

		for (a=0; a<Initial_I; a++)					//Set initial genotypes of cadavers
		{
			for (b=0; b<Initial_Strains; b++)
			{
				I[a][b]=0;
				//I[a][b]=gsl_rng_uniform_int(r, 2);
			/*	if (gsl_rng_uniform(r)<(double)b/(double)SNPs)
				{
					I[a][b]=0;
				}
				else
				{
					I[a][b]=1;
				}
			*/	//printf("I[%ld][%ld]=%lf\n", a, b, I[a][b]);
			}
			I[a][gsl_rng_uniform_int(r, Initial_Strains)]=1;			//Should this be Initial_Strains instead of SNPs?
		}



		printf("Seed=%d\n", atoi(argv[1]));
		printf("Years=%ld\n", Years);
		printf("phi=%lf\n", phi);
		printf("lambda=%lf\n", lambda);
		printf("CV=%lf\n", CV);
		printf("Initial_S=%ld\n", Initial_S);
		printf("Initial_I=%ld\n", Initial_I);
		printf("beta=%e\n", beta);
		printf("mu=%e\n", mu);
		printf("BaseReproduction=%e\n", BaseReproduction);
		printf("gamma=%e\n", gamma);
		printf("Model=%d\n", Model);


		for (j=0; j<Years; j++)
		{
			fflush(stdout);

			long int ExposureCounter=0;

			if (j!=0)
			{
				//printf("Total_S=%ld, TotalTotal_V=%ld\n", Total_S, TotalTotal_V);
//				if(Total_S==0 )//|| TotalTotal_V==0)
//				{
//					printf("Epidemic over\n\n");
//					TotalTotal_V=1;
//					//getchar();
//					break;
//				}
				//printf("TotalTotal_V=%ld\n", TotalTotal_V);

				//printf("Total_S=%ld\n", Total_S);

				Total_S=Initial_S;	//lambda*Total_S;
				Initial_S=0;

				//printf("Initial_I=%ld, gamma=%lf, ", Initial_I, gamma);

				TempOffspring=gsl_ran_binomial(r, gamma, Initial_I);
				Total_I=TempIndex+TempOffspring;

				//printf("TempOffspring=%ld\n", TempOffspring);

				for (a=TempIndex; a<Total_I; a++)
				{
					long int OldStrain=gsl_rng_uniform_int(r, Initial_I);

					double TempSum=0;
					for (b=0; b<Initial_Strains; b++)
					{

						I[a][b]=Old_I[OldStrain][b];

						TempSum+=I[a][b];
					}
					if (!DOUBLE_EQ(TempSum, 1))
					{
						printf("TempSum=%lf, j=%ld, a=%ld, OldStrain=%ld\n", TempSum, j, a, OldStrain);
					}
				}

				//if (j==6)
				//{
				//	for (b=0; b<Initial_Strains; b++)
				//	{
				//		if(j==7 ){printf("I[0][b]=%lf\n", I[0][b]);}
						//printf("I[0][%ld]=%lf\n", b, I[0][b]);
				//	}

				//}

				Initial_I=Total_I;

				//printf("Total_I=%ld\n", Total_I);

				CadaverCounter=0;

				ProbOfPredation=2*TypeThree_a*TypeThree_b*Total_S/(TypeThree_b*TypeThree_b + Total_S*Total_S);

				//printf("Predation=%lf\n", ProbOfPredation);


				beta=0;
				YearQuality=1;	//exp(gsl_ran_gaussian(r, 0.1));
				//YearQuality=exp(gsl_ran_gaussian(r, 0.1));
				for (a=0; a<Total_S; a++)
				{
					if (DOUBLE_EQ(S[a][0], Uninfected))
					{
						if (gsl_rng_uniform(r)<ProbOfPredation)			//Predation before reproduction
						{
							//printf("Here\n");
						}
						else
						{
							//printf("OrHere\n");
							TempOffspring=gsl_ran_poisson(r, YearQuality*(BaseReproduction+h[a]*lambda));
							Initial_S+=TempOffspring;
							beta+=(TempOffspring*h[a]);		//beta??
						}
					}
				}
				beta=beta/(double)Initial_S;
				//printf("beta=%e\n\n", beta);

				Initial_S++;				//Add in migration of a single host

				if(Initial_S>LargeNumber)
				{
					printf("Population size is out of hand\n");
					exit(EXIT_SUCCESS);
				}
			}

//			printf("Initial_S=%ld, Initial_I=%ld, Total_V=%ld, mean_transmission=%e, j=%ld\n", Initial_S, Initial_I, TotalTotal_V, beta, j);
			printf("%ld, %ld, %e, %ld\n", Initial_S, Initial_I, beta, j);
			TotalTotal_V=0;

			for (a=0; a<Initial_I; a++)
			{
				for (b=0; b<Initial_Strains; b++)
				{
					Old_I[a][b]=I[a][b];
				}
			}


			Total_S=Initial_S;
			Total_E=0;
			double Temp_h;

			//printf("Total_S=%ld\n", Total_S);

			double *E=malloc(Initial_S*sizeof(double));	if (E==NULL){printf("Too large\n");}


			double EndTime=1344;		//56 days * 24 hours/day;		//old value of endtime was double EndTime=2000;


			double Total_h=0;


			for (a=0; a<Initial_S; a++)
			{
			//printf("Here2\n");
				S[a][0]=-1;
	//printf("Here1\n");

				E[a]=1e100;
				h[a]=gsl_ran_gamma(r, 1/(CV*CV), beta*CV*CV);
				Total_h+=h[a];
			}

			///For speed, sort hetereogeneity values
			if (1)
			{
				//gsl_sort (h, 1, Initial_S);
				ShellSort(h, 0, Initial_S-1);
			}

			for (a=0; a<Initial_S; a++)
			{
				//printf("h[a]=%lf\n", h[a]);
			}


			double NextInfection;
			double TimeToNextEvent;
			double TypeOfNextEvent;

			double PredictedDeathTime;

			double ExposureRate;
			double RemovalRate;

			long int NextInfectedCat=0;


			double t=0;
			double WhichCat_Cutoff, WhichCat_Sum;		// values to determine which cat becomes exposed
			long int CaterpillarIndex;					// values to determine which cat becomes exposed

			long int StrainIndex;
			double ForceOfSecondaryInfection;
			long int WhichDeathTime;


			int NextInfectedChanged=1;
			double Temp;
			double ForceOfI;

			while(t<EndTime && (Total_E>0 || Total_I>0) && Total_S>0)				//while we're having an epidemic...
			{
				if (j==8)
				{
					for (a=0; a<Initial_Strains; a++)
					{
						//printf("S1[29][%ld]=%lf\n", a, S[29][a]);
					}
				}

				//if (Total_E%1000==0){printf("Total_E==%ld, ExposureCounter=%ld, t=%lf\n", Total_E, ExposureCounter, t);}

				////////////////Figure out what happens next//////////////////
				ExposureRate=Total_h*Total_I;
				RemovalRate=mu*Total_I;

				TimeToNextEvent=-log(gsl_rng_uniform(r))/(ExposureRate+RemovalRate);

				if ((NextInfectedChanged==1) && (t+TimeToNextEvent>0))
				{
					NextInfectedCat=gsl_stats_min_index(E, 1, ExposureCounter);
					//NextInfectedCat=gsl_stats_min_index(E, 1, Initial_S);
					NextInfectedChanged=0;
				}
				else if (NextInfectedChanged==2)
				{
					if (E[NextInfectedCat]>E[ExposureCounter-1])
					{
						NextInfectedCat=ExposureCounter-1;
					}
					NextInfectedChanged=0;
				}

				NextInfection=E[NextInfectedCat]-t;

				TypeOfNextEvent=gsl_rng_uniform(r)*(ExposureRate+RemovalRate);

				////////////////End figure out what happens next///////////////

				//if (j==6)
				//{
					//printf("TimeToNextEvent=%lf\n", TimeToNextEvent);
				//	printf("gsl_stats_min_index(E, 1, ExposureCounter)=%ld\n", gsl_stats_min_index(E, 1, ExposureCounter));
				//	printf("E[gsl_stats_min_index(E, 1, ExposureCounter)=%lf\n", E[gsl_stats_min_index(E, 1, ExposureCounter)]);
				//}

				//if (TimeToNextEvent>1e5 & NextInfection>TimeToNextEvent)
				//{
				//	break;
				//}

				if (NextInfection<TimeToNextEvent)			//Next event is a cat dying
				{
					NextInfectedChanged=1;
					//printf("1\n");
					Total_E--;
					Total_S--;
					Total_I++;




					//if(j==6)
					//{

					//	printf("E[NextInfectedCat]=%lf, NextInfectedCat=%ld, t=%lf\n", E[NextInfectedCat], NextInfectedCat, t);
					//	for (b=0; b<Initial_Strains; b++)
					//	{
					//		printf("S[NextInfectedCat][b]=%lf\n", S[NextInfectedCat][b]);
					//	}
					//	exit(EXIT_SUCCESS);
					//}

					for (a=0; a<Initial_Strains; a++)
					{
						I[Total_I-1][a]=S[NextInfectedCat][a];
						Total_V[CadaverCounter][a]=S[NextInfectedCat][a];
					//	if (j==8)
					//	{
					//		printf("S[%ld][%ld]=%lf\n", NextInfectedCat, a, S[NextInfectedCat][a]);
					//	}

					}
					CadaverCounter++;

					S[NextInfectedCat][0]=Dead;
					E[NextInfectedCat]=1e100;

					Total_h-=h[NextInfectedCat];
					h[NextInfectedCat]=0;


					t+=NextInfection;
				}

				else if (TypeOfNextEvent<ExposureRate)		//Next event is an exposure
				{
					//printf("2\n");

					StrainIndex=gsl_rng_uniform_int(r, Total_I);
					WhichDeathTime=gsl_rng_uniform_int(r, NumberOfDeaths);

					WhichCat_Cutoff=gsl_rng_uniform(r)*Total_h;
					WhichCat_Sum=0;
					for (a=0; a<Initial_S; a++)					//Figure out which caterpillar became infected
					{
						if (S[a][0]>=-2)					//Everything except dead caterpillars
						{
							WhichCat_Sum+=h[a];

							if (WhichCat_Sum>WhichCat_Cutoff)
							{
								CaterpillarIndex=a;
								break;
							}
						}
					}

					if (E[CaterpillarIndex]-t-TimeToNextEvent<100)
					{
						continue;
					}

					//PredictedDeathTime=t+TimeToNextEvent+DeathTime(I[StrainIndex], Initial_Strains, MeanKillTime, r);
					PredictedDeathTime=t+TimeToNextEvent+DeathHolder[WhichDeathTime];

					if (PredictedDeathTime>1e50)
					{
						NextInfectedChanged=0;
						continue;
					}


					NextInfectedChanged=2;			//recalculate time of next death

					if (DOUBLE_EQ(S[CaterpillarIndex][0], Uninfected))			//if the caterpillar was uninfected
					{
						//Swap caterpillar position to put the exposure time on the left of E
						for (a=0; a<Initial_Strains; a++)
						{
							S[CaterpillarIndex][a]=S[ExposureCounter][a];
						}
						E[CaterpillarIndex]=E[ExposureCounter];

						Temp_h=h[CaterpillarIndex];
						h[CaterpillarIndex]=h[ExposureCounter];

						Total_E++;

						///////////////////////////Begin to figure out which strains cause infection

						WithinHostDrift(TempHolder, I[StrainIndex], StrainHolder[WhichDeathTime], Initial_Strains, Model, r);
						for (a=0; a<Initial_Strains; a++)
						{
							S[ExposureCounter][a]=TempHolder[a];

							//if (ExposureCounter==24)
							{
								//printf("S[ExposureCounter][a]=%lf\n", S[ExposureCounter][a]);
							}
							//printf("TempHolder[%ld]=%lf\n", a, TempHolder[a]);
							//printf("I[%ld]=%lf\n", a, I[StrainIndex][a]);
							//printf("StrainHolder[%ld]=%lf\n", a, StrainHolder[WhichDeathTime][a]);
							//printf("S[%ld]=%lf\n", a, S[ExposureCounter][a]);
						}
						//for (a=0; a<SNPs; a++)
						//{
						//	S[ExposureCounter][a]=I[StrainIndex][a];		//Change this
						//}

						E[ExposureCounter]=PredictedDeathTime;

						h[ExposureCounter]=Temp_h;

						ExposureCounter++;
						TotalTotal_V++;

					}

					else			//S==Dead should never happen, so I shouldn't need to exclude it.
					{

						WithinHostDrift(TempHolder, I[StrainIndex], StrainHolder[WhichDeathTime], Initial_Strains, Model, r);

						double NewN=Params.N/(exp(Params.phi*(PredictedDeathTime-E[CaterpillarIndex])));
						ForceOfSecondaryInfection=NewN/(Params.N+NewN);

						//printf("NewDeathTime=%lf, OldDeathTime=%lf, ", PredictedDeathTime, E[CaterpillarIndex]);

						E[CaterpillarIndex]=E[CaterpillarIndex]+log(Params.N/(Params.N+NewN))/Params.phi;

						//printf("CombinedDeathTime=%lf\n", E[CaterpillarIndex]);


						for (a=0; a<Initial_Strains; a++)
						{
							S[CaterpillarIndex][a]=ForceOfSecondaryInfection*TempHolder[a]+(1-ForceOfSecondaryInfection)*S[CaterpillarIndex][a];
						}



///ErrorChecking
						double TempSum=0;
						for (b=0; b<Initial_Strains; b++)
						{
							TempSum+=S[CaterpillarIndex][b];
						}
						if (!DOUBLE_EQ(TempSum, 1))
						{
							printf("TempSum=%lf, j=%ld, CaterpillarIndex=%ld\n", TempSum, j, CaterpillarIndex);
							exit(EXIT_SUCCESS);
						}
/////ErrorChecking



						//printf("Coinfection: Before S[%ld]=%lf\n", CaterpillarIndex, S[CaterpillarIndex]);
/*						if (E[CaterpillarIndex]<PredictedDeathTime)		//original infection is going to dominate
						{
							ForceOfSecondaryInfection=(double)Params.N/exp(Params.phi*(PredictedDeathTime-E[CaterpillarIndex]));
							ForceOfSecondaryInfection=ForceOfSecondaryInfection/(ForceOfSecondaryInfection+Params.N);
							for(a=0; a<Initial_Strains; a++)
							{
								S[CaterpillarIndex][a]=ForceOfSecondaryInfection*TempHolder[a]+(1-ForceOfSecondaryInfection)*S[CaterpillarIndex][a];
								//S[CaterpillarIndex][a]=ForceOfSecondaryInfection*I[StrainIndex][a]/(ForceOfSecondaryInfection+Params.N)+Params.N*S[CaterpillarIndex][a]/(ForceOfSecondaryInfection+Params.N);
							}
						}
						else
						{
							ForceOfSecondaryInfection=(double)Params.N/exp(Params.phi*(E[CaterpillarIndex]-PredictedDeathTime));
							ForceOfSecondaryInfection=ForceOfSecondaryInfection/(ForceOfSecondaryInfection+Params.N);
							//ForceOfSecondaryInfection=Params.N*exp(-Params.phi)/(E[CaterpillarIndex]-Params.N);
							for (a=0; a<Initial_Strains; a++)
							{
								S[CaterpillarIndex][a]=ForceOfSecondaryInfection*S[CaterpillarIndex][a]+(1-ForceOfSecondaryInfection)*TempHolder[a];
								//S[CaterpillarIndex][a]=ForceOfSecondaryInfection*S[CaterpillarIndex][a]/(ForceOfSecondaryInfection+Params.N)+Params.N*I[StrainIndex][a]/(ForceOfSecondaryInfection+Params.N);
							}
							E[CaterpillarIndex]=PredictedDeathTime;
						}
*/						//printf("Coinfection: After S[%ld]=%lf\n", CaterpillarIndex, S[CaterpillarIndex]);
						//getchar();
					}

					t+=TimeToNextEvent;
				}

				else
				{
					NextInfectedChanged=0;

					StrainIndex=gsl_rng_uniform_int(r, Total_I);

					for (a=0; a<Initial_Strains; a++)
					{
						I[StrainIndex][a]=I[Total_I-1][a];
					}

					Total_I--;
					t+=TimeToNextEvent;
				}
			}

			//printf("%ld %ld %ld %lf\n", Total_S, Total_I, Total_E, t);

/*			long int NumberOfStrains=0;
			for (a=0; a<MaxStrains; a++)
			{
				if (Total_V[a]>0)
				{
					NumberOfStrains++;
					//printf("V[%ld]=%ld\n", a, Total_V[a]);
				}
			}
*/			//printf("NumberOfStrains=%ld\n", NumberOfStrains);
			//printf("Survivors=%ld\n", Total_S);
			//printf("%ld\n", Total_S);
	///////////////////Calculate pi/////////////////////////////
/*
			double pi=0;
			TotalTotal_V=0;
			for(a=0; a<MaxStrains; a++)
			{
				TotalTotal_V+=Total_V[a];
			}
			for (a=0; a<MaxStrains; a++)
			{
				for (b=0; b<a; b++)
				{
					if ((Total_V[a]>0) && (Total_V[b]>0))
					{
						pi+=(double)(Differences(&Genotype[a][0], &Genotype[b][0], SNPs)*(double)Total_V[a]*(double)Total_V[b])/(double)((TotalTotal_V*TotalTotal_V*SNPs));
					}
				}
			}
*/

			for (a=0; a<Initial_S; a++)  //To account for the caterpillars that are still exposed when the epidemic ends
			{
				if (S[a][0]>Uninfected+.001)
				{
					for (b=0; b<Initial_Strains; b++)
					{
						Total_V[CadaverCounter][b]=S[a][b];
					}
					CadaverCounter++;
				}
			}


			/////////////Begin find Allele Frequencies


			if(j==Years-1)
			{
				if ((CadaverCounter<2000) || (Initial_S<10000))	//make these better
				{
					Years++;
				}
				else
				{
				}
			}






			if (j==Years-1)
			{

				long int c;

				long int RandomSample[2000];
				long int Temp[CadaverCounter];

				for (i=0; i<CadaverCounter; i++)
				{
					Temp[i]=i;
				}

				gsl_ran_choose (r, RandomSample, 2000, Temp, CadaverCounter, sizeof(long int));


//				for (c=0; c<CadaverCounter; c++)
				for (c=0; c<2000; c++)
				{
					for (b=0; b<Initial_Strains; b++)
					{
						printf("%.3f ", Total_V[RandomSample[c]][b]);
						//SNP_Freq[c][a]+=Total_V[c][b]*TempGenotype[b][a];
						//printf("Increment=%lf\n", Total_V[c][b]*TempGenotype[b][a]);
					}
					printf("\n");
				}
			}







			TempIndex=0;
			//for (b=0; b<Initial_Strains; b++)
			//printf("Total_V_or_I[0][0]=%lf\n", I[0][b]);
			for (a=0; a<CadaverCounter; a++)
			{

				TempOffspring=gsl_ran_poisson(r,phi);


				for (c=0; c<TempOffspring; c++)
				{
					for (b=0; b<Initial_Strains; b++)
					{
						I[TempIndex][b]=Total_V[a][b];
					}
					TempIndex++;
				}
			}

			TempStrain=gsl_rng_uniform_int(r, Initial_Strains);
			for (a=0; a<Initial_Strains; a++)
			{
				if (a==TempStrain)
				{
					I[TempIndex][a]=1;
				}
				else
				{
					I[TempIndex][a]=0;
				}
			}
			TempIndex++;
			CadaverCounter++;


			free(E);

		}
	}

	//printf("%ld\n", time(NULL)-StartTime);

	free(h);
	free(S);
	free(Total_V);
	exit(EXIT_SUCCESS);

}
