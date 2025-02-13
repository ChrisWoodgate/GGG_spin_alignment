/*==============================================/
/                                               /						
/ spin_functions header file                    /
/ Original code created by C Woodgate - 2018    /
/                                               /
/==============================================*/

/* Function prototypes */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include "mt19937ar.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void B_Field(double B[], double theta, double phi);

double BFun(double R, double alpha);

double CFun(double R, double alpha);

double complex Hopping(double q_x, double q_y, double q_z, double del_x, double del_y, double del_z);

double realdot(double *vec1, double *pos, int j);

int pos_idx(int i, int j);

void k_eff_field(double *EffFieldK, double *Dir, double *position, int N, double alpha, double Dipol, int k_max);

void R_eff_field(double *EffFieldR, double *Dir, double *position, int N, double alpha, 
				double Dipol, int L_max, double J1, int NCells);

void Eff_field(double BB[3], double *EffField, double *EffFieldK, double *EffFieldR, double *Dir, int N, double alpha, double Dipol, double Mag, double Zeeman);

double Energy(double BB[3], double *EffField, double *Dir, double *energies, 
			int N, double Mag, double Zeeman);

void cart_to_polar(double *cart, double polar[]);

void polar_to_cart(double *cart, double polar[]);

void Noise(double zhat[], double size);

void Peturb(double *Dir, int N, double size);
