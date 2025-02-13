/*==============================================/
/                                               /						
/ Data handling .h file                         /
/ Original code by C Woodgate - 2018            /
/												/
/==============================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>

void read_input(int *N_ions, double *alpha, int *k_max, int *L_max, double *theta, double *phi); 

void read_positions(long N_ions, double *positions, int NCells);

int read_args(int argc, char **argv, int *spin_movement, int *pertubation, int *NSteps, int *NCells, double *field_strength);
