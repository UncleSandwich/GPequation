/************************************************************/
/*  This code is aimmed to calculate the energies of the ground and excied states of Gross-Pitaevskii equation   */
/*  Euqation is originated from the BEC many-body problems   */
/*  Using the standard imaginary time propagation method     */
/*  The original code is written by W. Zang, on Oct. 15, 2003  */
/*  Code rewritten by Hangxi, Li */
/*  Add calculation of the nth state's energy with the modification of wave_func() and evolution()*/
/************************************************************/

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include "engine.h"

/* Define the globel constants */
#define PI 3.14159265    /* circumference ratio, a common constant */
#define NZ 625         /* the total number of points along z-direction */
#define ZLEN 32.0       /* the range (length) of function */

#define LIMIT 1.5e-7   /* the threshold value for comparison  */
#define PRINT 100000     /* the number of rounds to output the function values  */
#define BREAK 50000000    /* used to prevent the dead loop */

double dt = 0.0000001;  /* imaginary time iteration dt may be equal to Dz*Dz */
const double t = 0.0000001;

const double Dz = ZLEN / (NZ - 1.);  /* the spacing of points */
const double Zmid = (NZ - 1.) / 2;   /* the middle point */
/* ------------------------------------------------------------------------*/

/* Declare global variables */
double energy;   /* the state energy */
double pb_r[NZ], pb_i[NZ]; /* the real part and the imaginary part of wave function */
double pb[NZ]; // The probability of wave function
double lamda;    /* the parameter that used to find the excited states */
double a[2],b[2];//,c[2]; // Three random arrays
int w; //w=0 or 1 to determine calculating even or odd

//The variables shall be controled
double g;   /* a constant, g = 4*pi*hbar^2*a_s/m , in BEC problem */
double gini = 10.0; // the initial value of g
double gup = 10.0; // the upbound value of g
double ginter = 1.0; // the interval of the loop of g

double ini_val = 9.0;   /* the initial value of lamda */
double upbound = 9.1;  /* the upbound value of lamda */
double interval = 0.01;    /* the interval of the loop of lamda */


/* Declare functions */
void random_gener();   // generate random number for initial wave function
void ini_wave();    /* function to create the initial wave function, Guassian form */
void evolution();   /* function to evolve the wave funcion with time iteration*/
void evolution_origi(); /* function to evolve the wave funcion with original Hamiltonian operator*/
void energy_cal();  /* function to calculate the ground state energy of present wave function */
void normalize();   /* function to normalize the wave functin after each time's evolution */

void wave_func();  /* function to evolve, normalize the wave function and calculate the energies by calling existed functions,
				   and compare the wave function evolved with the previous one before a time iteration to decide whether to save the result */

void save_func();    /* function to save the final wave function */
void save_posi();    /* function to save the position of points */
void save_result();  /* function to save the energies*/
/* ------------------------------------------------------------------------*/

int main()
{
	save_posi();     /* save the position of points */
	//ini_wave();      /* create an guassian wave funcion as the initial */
	//energy_cal();    /* calculate the ground state energy of this initial wave function*/
	/*save_energy();   /* save the initial energy */
	
	
        for(g = gini; g <= gup; g += ginter)
	{
	
		 wave_func();
	}	

	system("pause");
	return 0;
}
/*--------------------------------------------------------------------------*/

/* function to save the position of points */
void save_posi()
{
	char file[30]; /* used for naming the files */
	FILE *fp;

	/* name the file storing the position of points */
	strcpy(file, "z");
	strcat(file, ".dat");
	/* input the position of points along z-direction */
	fp = fopen(file, "w");
	int i;
	for (i = 0; i<NZ; i++)
		fprintf(fp, "%f\n", (i - Zmid)*Dz);
	fclose(fp);
}
/*--------------------------------------------------------------------------*/

// generate random number for initial wave function
void random_gener()
{
	int i;
	
	srand((unsigned)time(NULL));
	
	do{
		for (i = 0; i<2; i++)
        	{
                	a[i] = (double)(rand() % 300 + 5);
              		//  c[i] = (double)(rand() % 200);
                	b[i] = (double)(rand() % 300 + 5);
        	}
	}while( a[1] == 0 || a[0] == 0 || b[0] == 0 || b[1] == 0 );

	//for (i = 0; i<2; i++)
	//{	
	//	printf("%f %f \n", a[i], b[i]);
	//}
		
}

/* function to create the initial wave function, Guassian form */
void ini_wave()
{
	double C = pow((1. / PI), 0.25); /* define the constant coefficient in guassian function*/
	double z;
	int i;

	//srand((unsigned)time(NULL));
	
	printf("%f %f\n", a[0], b[0]);
	/*for (i = 0; i < NZ; i++)
	{
		z = (i - Zmid)*Dz;
		pb_r[i] = C*exp(-0.5*z*z / 10);
		pb_i[i] = C*exp(-0.5*z*z / 10);
	}*/
	if (w == 0)
	{
		for (i = 0; i<NZ; i++)
		{	
			z = (i - Zmid)*Dz;
			pb_i[i] = cos(z/b[0]);
			pb_r[i] = cos(z/b[1]);
		}
	}
	else if (w == 1)
	{
		for (i = 0; i<NZ; i++)
		{
			z = (i - Zmid)*Dz;
			pb_i[i] = sin(z/a[0]);
			pb_r[i] = sin(z/a[1]);
		}
	}
	else 
	{
		for (i = 0; i<NZ; i++)
		{
			pb_i[i] = rand() % (int)(a[0]) + 1;
			pb_r[i] = rand() % (int)(b[1]) + 1;
		}
	}

	normalize();
	/*for (i = 0; i<NZ; i++)
		printf("%f ", pb_r[i]);
	printf("\n");
	for (i=0;i<NZ;i++)
		printf("%f ", pb_i[i]);
	printf("\n");
	system("pause");*/	
	/*for (i = 0; i < NZ; i++)
	{
		pb[i] = pb_r[i] * pb_r[i] + pb_i[i] * pb_i[i];
	}*/

	//char file[60]; /* used for naming the files */
	//FILE *fp;
	///* name the file storing the initial wave function */
	//strcpy(file, "initial wavefunction");
	//strcat(file, parameters);
	//strcat(file, ".dat");
	///* input the values of wave function */
	//fp = fopen(file, "w");
	//for (i = 0; i < NZ; i:++)
	//	fprintf(fp, "%25.15e\n", pb[i]);
	//fclose(fp);
}
/*--------------------------------------------------------------------------*/

/* function to calculate the ground state energy of present wave function */
void energy_cal()
{
	double z;
	double V; /* the harmonic potential */
	double n; /* the square amplitude of wave function at a point */
	/* discrete the wave function, define the value of the nearest, and the second nearest points around a point */
	double pc, pl, pr, pll, prr; /* c-centre, l-left, r-right, ll-second left, rr-second right */
	double ekz_r, ekz_i; /* the real and imaginal part of the double differential of the wave function, for calculating kinetic energy */

	energy = 0.0; /* initialize the variable */

	int i;
	#pragma omp for
	for (i = 0; i<NZ; i++)
	{
		z = (i - Zmid)*Dz;
		V = 0.5*z*z;
		n = pb_r[i] * pb_r[i] + pb_i[i] * pb_i[i];

		/* the real part of wave function of those five points */
		pc = pb_r[i];
		pll = (i <= 1)      ? 0. : pb_r[i - 2];
		prr = (i >= NZ - 2) ? 0. : pb_r[i + 2];
		pl = (i == 0)       ? 0. : pb_r[i - 1];
		pr = (i == NZ - 1)  ? 0. : pb_r[i + 1];
		ekz_r = (8.*(pr - pl) - (prr - pll)) / 12. / Dz; /* the real part of double differential */

		/* the imaginary part of wave function of those five points */
		pc = pb_i[i];
		pll = (i <= 1)      ? 0. : pb_i[i - 2];
		prr = (i >= NZ - 2) ? 0. : pb_i[i + 2];
		pl = (i == 0)       ? 0. : pb_i[i - 1];
		pr = (i == NZ - 1)  ? 0. : pb_i[i + 1];
		ekz_i = (8.*(pr - pl) - (prr - pll)) / 12. / Dz; /* the imaginal part of double differential */

		/* energy without a factor Dz */
		energy += 0.5*(ekz_r*ekz_r + ekz_i*ekz_i) + V*n + 0.5*g*n*n; /* mass, hbar in kinetic term have been ignored */
	}
	energy *= Dz;
}
/*--------------------------------------------------------------------------*/

/* function to save the energies*/
void save_result()
{
	char file[60]; /* used for naming the files */
	FILE *fp;
	
	strcpy(file, "Energy and lamda and g");
	strcat(file, ".dat");
	/*input the value of energies */
	fp = fopen(file, "a");
	fprintf(fp, "%f   %12.11f  %f \n", lamda, energy, g);
	fclose(fp);
	
}
/*--------------------------------------------------------------------------*/

/* function to evolve, normalize the wave function and calculate the energies by calling existed functions,
and compare the wave function evolved with the previous one before a time iteration to decide whether to save the result */


void wave_func()
{
	double tem_energy; /* a variable used to record the temporary value of previous energy to make a comparison */
	double diff_energy; /* used to record the difference of energy between evolved wave function and the previous one */
	double diff_wavefunc = 0.; /* used to record the cumulated difference of values between evolved wave function and the previous one */
	double tem_pb_r[NZ] = { 0. }, tem_pb_i[NZ] = { 0. }; /* used to record the temporary values of previous wave function to make a comparison */
	double tem_pb[NZ]; //The probability of pre-evolved wave functions
	int round; /* record the number of rounds */

	int i, j;
	for (lamda = ini_val; lamda < upbound; lamda += interval)
	{
		round = 1;
		dt = 0.0000001;  /* imaginary time iteration dt may be equal to Dz*Dz */
		j = 0;
		w = 0;
		w = (int) lamda % 2; 

		random_gener();

		ini_wave();
		do{
			tem_energy = energy; /* record the energy */
			/* record the wave function */
			for (i = 0; i<NZ; i++)
			{
				tem_pb_r[i] = pb_r[i];
				tem_pb_i[i] = pb_i[i];
			}

			//calculate the probability of pre-evolved wave function
			for (i = 0; i < NZ; i++)
			{
				tem_pb[i] = tem_pb_r[i] * tem_pb_r[i] + tem_pb_i[i] * tem_pb_i[i];
			}

			/* get the wave function and energy after one iteration evolution */
			evolution();
			normalize();
			energy_cal();

			/* calculate the difference between two energies */
			diff_energy = fabs(energy - tem_energy);

			//calculate the probability of evolved wave function
			for (i = 0; i < NZ; i++)
			{
				pb[i] = pb_r[i] * pb_r[i] + pb_i[i] * pb_i[i];
			}

			/* calculate the cumulated difference between two wave functions */
			diff_wavefunc = 0.0; /* initialize the variable */
			for (i = 0; i < NZ; i++)
				diff_wavefunc += fabs(tem_pb[i] - pb[i]);

			/* add a little of pre-evolved wave function to the post-evolved wave function to prevent the chaotic phenomena, edited on June 1st, 2015 by Li Hangxi*/
			/*double a = 0.01;
			for (i = 0; i < NZ; i++)
			{
				pb_r[i] += a*tem_pb_r[i];
				pb_i[i] += a*tem_pb_i[i];
			}

			/* decide whether to output the temporary results according to the numbeer of rounds */
			if (round % PRINT == 0)
			{
				if (dt <= t*500)
					dt += 0.1*t; //As the evolution proceeds, the iteration time is increasing
				printf("round=%5d	%12.11f	%e	%e  %f\n", round, energy, diff_energy, diff_wavefunc, lamda);
			}
			/* decide whether to break the loop */
			if (round % BREAK == 0)
			{
				printf("Break\n");
				break;
			}
				
			if ( fabs(diff_wavefunc - LIMIT) < 1.0e-9)
				j++;
			round++;
		} while (j != 2000 && diff_wavefunc > LIMIT);

		
		/* save the results: energy and wave function */
		save_result();
		save_func();

		/* save the results: temporary results */
		char file[60]; 
		FILE *fp;		
		strcpy(file, "result.txt");
		fp = fopen(file, "a");
		fprintf(fp, "NZ=%d, ZLEN=%f, LIMIT=%d, PRINT=%d, BREAK=%d, \n", NZ, ZLEN, LIMIT, PRINT, BREAK);
		fprintf(fp, "round=%d, j=%d, t=%e, dt=%e,\n", round,j,t,dt);
		fprintf(fp, " a[0]=%f  a[1]=%f  b[0]=%f  b[1]=%f,  lamda=%f,  g=%f \n", a[0], a[1], b[0], b[1], lamda, g);
		fprintf(fp, "energy=%e, diff_energy=%e, diff_wavefunc=%e\n\n", energy, diff_energy, diff_wavefunc);
		fclose(fp);
	}


}
/*--------------------------------------------------------------------------*/


/* direct function to evolve the wave funcion with time iteration, written at June 12, 2015 by Li Hangxi*/
void evolution()
{
	double z;
	double tem_pb_r[NZ], tem_pb_i[NZ]; /* used to record wave function before any operation */
	double pc, pl, pr, pll, prr; /* c-centre, l-left, r-right, ll-second left, rr-second right */
	double pci, pli, pri, plli, prri; // the imaginary part of wave function with the same meanings of subscript

	int i;
	
	for (i = 0; i<NZ; i++)
	{
		tem_pb_r[i] = pb_r[i];
		tem_pb_i[i] = pb_i[i];
	}

	/* the time evolution of the wavefunction */
	#pragma omp for
	for (i = 0; i<NZ; i++)
	{
		z = (i - Zmid)*Dz;

		/* the real part of wave function */
		pc = pb_r[i];
		pl  = (i == 0)       ? 0.0 : tem_pb_r[i - 1];
		pr = (i == NZ - 1)   ? 0.0 : tem_pb_r[i + 1];
		pll = (i <= 1)       ? 0.0 : tem_pb_r[i - 2];
		prr = (i >= NZ - 2)  ? 0.0 : tem_pb_r[i + 2];

		/* the imaginary part of wave function */
		pci  = pb_i[i];
		pli  = (i == 0)        ? 0.0 : tem_pb_i[i - 1];
		pri  = (i == NZ - 1)   ? 0.0 : tem_pb_i[i + 1];
		plli = (i <= 1)        ? 0.0 : tem_pb_i[i - 2];
		prri = (i >= NZ - 2)   ? 0.0 : tem_pb_i[i + 2];

		/* The real part of evolved wave function*/
		pb_r[i] = pc + (-dt) * 
			(
			0.25 / (Dz*Dz*Dz*Dz)*(prr - 4 * pr + 6 * pc - 4 * pl + pll) 
			+ 0.5 / (Dz*Dz)*(pr - 2 * pc + pl)*(-z*z - 3 * g*(pc*pc + pci*pci) + 2 * lamda)
			- 0.5*g / Dz*(pc*(pr*pr - pri*pri + pl*pl - pli*pli - 2 * pr*pl + 2 * pri*pli) + 2 * pci*(pr*pri + pl*pli - pri*pl - pr*pli))
			- 0.5*g / (Dz*Dz)*pc*((pr*pr + pri*pri) + (pl*pl + pli*pli) - 2 * (pl*pr + pli*pri))
			- 0.5*g / (Dz*Dz)*((pc*pc - pci*pci)*(pr - 2 * pc + pl) - 2 * pc*pci*(-pri + 2 * pci - pli))
			- 0.5*z / Dz*(pr - pl)
			+ pc*(-0.5 + lamda*lamda + 0.25*z*z*z*z - lamda*z*z + (g*z*z - 2 * lamda*g)*(pc*pc + pci*pci) + g*g*(pc*pc + pci*pci)*(pc*pc + pci*pci))
			);

		/* The imaginary part of evolved wave function*/
		pb_i[i] = pci + (-dt)*
			(
			0.25 / (Dz*Dz*Dz*Dz)*(prri - 4 * pri + 6 * pci - 4 * pli + plli) 
			+ 0.5 / (Dz*Dz)*(pri - 2 * pci + pli)*(-z*z - 3 * g*(pc*pc + pci*pci) + 2 * lamda)
			- 0.5*g / Dz*(2 * pc*(pr*pri + pl*pli - pri*pl - pr*pli) - pci*(pr*pr - pri*pri + pl*pl - pli*pli - 2 * pr*pl + 2 * pri*pli))
			- 0.5*g / (Dz*Dz)*pci*((pr*pr + pri*pri) + (pl*pl + pli*pli) - 2 * (pl*pr + pli*pri))
			- 0.5*g / (Dz*Dz)*(2 * pc*pci*(pr - 2 * pc + pl) - (pc*pc - pci*pci)*(-pri + 2 * pci - pli))
			- 0.5*z / Dz*(pri - pli)
			+ pci*(-0.5 + lamda*lamda + 0.25*z*z*z*z - lamda*z*z + (g*z*z - 2 * lamda*g)*(pc*pc + pci*pci) + g*g*(pc*pc + pci*pci)*(pc*pc + pci*pci))
			);

	}

}
/*--------------------------------------------------------------------------*/

/* function to normalize the wave functin after each time's evolution */
void normalize()
{
	double norm = 0.0; /* the total value of wave square amplitude */
	double factor = 0.0; /* the factor used to multiply to make wave function normalized*/

	int i;
	for (i = 0; i<NZ; i++)
		norm += pb_r[i] * pb_r[i] + pb_i[i] * pb_i[i];

	factor = sqrt(1. / (norm*Dz));
	for (i = 0; i<NZ; i++)
	{
		pb_r[i] *= factor;
		pb_i[i] *= factor;
	}
}
/*--------------------------------------------------------------------------*/

/* function to save the final wave function */
void save_func()
{
	char file[60]; /* used for naming the files */
	FILE *fp, *fpr, *fpi;
	int i;
	char string[20],parameters[7];
	snprintf(string, 20, "%5.3f", lamda);
	snprintf(parameters, 6, "%5.3f", g);

	/* output the wavefunction */
	strcpy(file, "wave lamda=");
	strcat(file, string);
	strcat(file, "g=");
	strcat(file, parameters);
	strcat(file, ".dat");
	/* input the values of wave function */
	fp = fopen(file, "w");
	//fprintf(fp, "lamda=%12.5f\n", lamda);
	for (i = 0; i < NZ; i++)
		fprintf(fp, "%25.15e\n", pb[i]);
	fclose(fp);

	/* output the real part of wavefunction */
	strcpy(file, "real lamda=");
	strcat(file, string);
	strcat(file, "g=");
	strcat(file, parameters);
	strcat(file, ".dat");
	fpr = fopen(file, "w");
	for (i = 0; i < NZ; i++)
		fprintf(fpr, "%25.15e\n", pb_r[i]);
	fclose(fpr);

	/* output the imaginary part of wavefunction */
	strcpy(file, "imaginary lamda=");
	strcat(file, string);
	strcat(file, "g=");
	strcat(file, parameters);
	strcat(file, ".dat");
	fpi = fopen(file, "w");
	for (i = 0; i < NZ; i++)
		fprintf(fpi, "%25.15e\n", pb_i[i]);
	fclose(fpi);

}
