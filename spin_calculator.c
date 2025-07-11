/*================================================/
/                                                 /				
/ Code for calculating spin canting in GGG        /
/ Original code created by C. D. Woodgate - 2018  /
/                                                 /
/================================================*/

#include "data.h"
#include "spin_functions.h"

int main(int argc, char **argv)
{
	/* Parameters to be read in from input file */
	int N_ions, N_unit; // Number of ions in original cell and in the system
	double alpha_ew;    // Value of alpha parameter to be used in Ewald summation
	int L_max;          // L_max and 
	int k_max;          // k_max are also parameters used in the Ewald summation
	double Theta;       // Theta
	double Phi;         // and Phi set direction of external field
	double Mag;         // External Field (T)

	/* Parameters to be read in from command line and default values  */
	int spin_movement      =    0; // If 1, Spins moved one by one. If 0, moved together.
	int pertubation        =    0; // If 1, Pertubation turned on. If 0, turned off.
	int N_iter             =   10; // N_iter
	double external_field  = 2400; // External Field (mT)
	int NCells             =    1; // Number of unit cells in x, y, z directions
	
	/* Read in from input file */
	read_input(&N_unit ,&alpha_ew, &k_max, &L_max, &Theta, &Phi);

	read_args(argc, argv, &spin_movement, &pertubation, &N_iter, &NCells, &external_field);

	N_ions = NCells*NCells*NCells*N_unit;
	Mag = 0.001*external_field;

	/* Allocate memory for positions array and read in from file */
	double *positions;
	positions = (double *)malloc(N_ions*3*sizeof(double));

	read_positions(N_ions, positions, NCells);
	
	/* Define relevant physical constants to be used in calculation */
	double NNex = 0.107 * 0.25 * 7.0;
	double Dipol = 0.0453 * 0.25 * 7.0;
//	double DddSpin = 0.0453 * pow(6, 1.5) * 0.25 * 7.0;
//	double NN1ex = Dipol;
//	double eVconvert = 1.38e-23 / 1.6e-22;
//	double lattspace = 12.349;
//	double InvAngst = 8.0/lattspace;
	double Zeeman = 9.27400968e-24 / 1.38e-23;
	double J1 = NNex;

	/* Allocate memory for storing contributions to field at each timestep */
	double *EffFieldK, *EffFieldR, *Dir, *EffFieldTot, *dS, *energies;
	Dir          = (double *)malloc(N_ions*3*sizeof(double));
	EffFieldK    = (double *)malloc(N_ions*3*sizeof(double));
	EffFieldR    = (double *)malloc(N_ions*3*sizeof(double));
	EffFieldTot  = (double *)malloc(N_ions*3*sizeof(double));
	dS           = (double *)malloc(N_ions*3*sizeof(double));
	energies     = (double *)malloc(N_ions*sizeof(double));

	/* Set up variable to store energy and make B-Field */
	double E_step;
	double B[3];
	B_Field(B, Theta, Phi);

	/* Output information to screen */
	printf("B Field = {%lf, %lf, %lf}, Mag = %lf, kmax = %d, Lmax = %d, alpha = %lf , N_ions = %d \n", 
			B[0], B[1], B[2], Mag, k_max, L_max, alpha_ew, N_ions);
	
	/* Align spins with external field to begin with */
	for(int i=0; i<N_ions; i++)
	{
		for(int j=0; j<3; j++)
		{
			Dir[pos_idx(i,j)] = B[j];
		}
	}

    int seed = time(NULL);
  	init_genrand(seed);
	
	/*=======================================/
	/    Outer loop over iterations          /
	/=======================================*/
	for(int i_step=0; i_step<N_iter; i_step++)
	{
		if(pertubation == 1 && i_step<5)
		{
			Perturb(Dir, N_ions, 0.1);
		}

		/*=======================================/
		/       Loop over ions in system         /
		/=======================================*/
		int ion;
		#pragma omp parallel for default(shared) private(ion)
		for(ion=0; ion<N_ions; ion++)
		{
			/* Do k-space calculation */
			local_k_eff_field(EffFieldK, Dir, positions, N_ions, alpha_ew, Dipol, k_max, ion);
			/* Do real-space calculation */
			local_R_eff_field(EffFieldR, Dir, positions, N_ions, alpha_ew, Dipol, L_max, J1, NCells, ion);
			/* Combine the two */
			local_Eff_field(B, EffFieldTot, EffFieldK, EffFieldR, Dir, N_ions, alpha_ew, Dipol, Mag, Zeeman, ion);

			if (spin_movement ==1)
			{
				//Move spins one by one
				Dir[  3*ion  ] = -EffFieldTot[  3*ion  ];
				Dir[3*ion + 1] = -EffFieldTot[3*ion + 1];
				Dir[3*ion + 2] = -EffFieldTot[3*ion + 2];
				
				// Normalise field
				normalise(Dir, ion);
			}
		} //Ion loop

		/* Compute Energies for each ion */
		E_step = Energy(B, EffFieldTot, Dir, energies, N_ions, Mag, Zeeman);
		
		/* File for writing intermediate steps to */	
		char filename[sizeof "outputs/configuration_cells01_mag0.00_steps000.csv"];
		sprintf(filename,"outputs/configuration_cells%02d_mag%3.2lf_steps%03d.csv",N_ions,Mag, i_step+1);
		FILE *outfile;
		outfile = fopen(filename,"w");
		if (outfile == NULL)
		{
			printf("Error opening output file\n");
			exit(EXIT_FAILURE);
		}

		/*=======================================/
		/       Move spins                       /
		/=======================================*/
		for(int i=0; i<N_ions; i++)
		{
			if (spin_movement == 0)
			{	
				//Move spins one by one
				Dir[  3*i  ] = -EffFieldTot[  3*i  ];
				Dir[3*i + 1] = -EffFieldTot[3*i + 1];
				Dir[3*i + 2] = -EffFieldTot[3*i + 2];
				
				// Normalise field
				normalise(Dir, i);
			}
			
			double mod_dS2 = 0;
			double mod_dS;	

			for(int j=0; j<3; j++)
			{
				int idx = pos_idx(i, j);
				/* New Field values */
				/* dS values */
				dS[idx] = Dir[idx] - B[j];
				mod_dS2 += dS[idx]*dS[idx]; 				
			}
			mod_dS = sqrt(mod_dS2);			

			/* Write to intermediate file */
			fprintf(outfile, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n", Dir[3*i], Dir[3*i+1], 
					Dir[3*i+2], dS[3*i], dS[3*i+1], dS[3*i+2], mod_dS, energies[i]);
				
			if(i==0)
			{

				double rotation[3];
				Cross(&Dir[3*i], B, rotation);
				normalise(rotation, 0);	
				printf("S_0 = (%lf, %lf, %lf), Rotation about (%5.2lf, %5.2lf, %5.2lf), Energy = %lf \n", Dir[3*i], Dir[3*i+1], Dir[3*i +2], rotation[0], rotation[1], rotation[2], E_step);
			}
		} //Ion loop
		fclose(outfile);
	} /* Iterations */
	
	/*=======================================/
	/    Clean up                            /
	/=======================================*/
	free(positions);
	free(Dir);
	free(EffFieldK);
	free(EffFieldR);
	free(EffFieldTot);
	free(dS);
	free(energies);
	return 0;
}
