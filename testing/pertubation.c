#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spin_functions.h"

int main(void)
{
	int N = 4;
	int M = 2;
	int i, j;
	double *spin;
	spin = (double *)malloc(3*sizeof(double));
	for(i=0; i<N; i++)
	{
		for(j=0; j<M; j++)
		{
			double theta = ((float)j)*M_PI/2 + M_PI/4;
			double phi = ((float)i)*M_PI/2;
			double polar[2] = {theta, phi};
			polar_to_cart(spin, polar);
			printf("Spin is (%lf, %lf, %lf)\n", spin[0], spin[1], spin[2]);
			Peturb(spin, 1, 0.01);
			printf("After pertubation, spin is (%lf, %lf, %lf)\n", spin[0], spin[1], spin[2]);
		}
	}

	free(spin);
	return EXIT_SUCCESS;
}
