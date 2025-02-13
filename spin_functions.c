/*==============================================/
/												/						
/ spin_functions module							/
/ Original code created by C Woodgate - 2018	/
/												/
/==============================================*/

#include "spin_functions.h"


void B_Field(double B[], double theta, double phi)
{
	/*Function to return a vector for given theta, phi*/
	double sin_theta = sin(theta);
	B[0] = sin_theta * cos(phi);
	B[1] = sin_theta * sin(phi);
	B[2] = cos(theta);
	return;
}

double BFun(double R, double alpha)
{
	/*B function as defined in Ewald summation reference*/
	double Rm3 = 1.0/(R * R * R);
	double output = Rm3*erfc(alpha*R);
	output += Rm3*(2.0 * alpha * R)/(sqrt(M_PI))*exp(-alpha*alpha*R*R);
	return output;
}

double CFun(double R, double alpha)
{
	/*C function as defined in Ewald summation reference*/
	double Rm5 = pow(R, -5.0);
	double output = Rm5*(3.0 * erfc(alpha * R));
	output += Rm5*((2.0 * alpha * R)/(sqrt(M_PI))*(3.0 + 2.0* alpha*alpha * R*R) * exp( - alpha*alpha * R*R));
	return output;
}

double complex Hopping(double q_x, double q_y, double q_z, double del_x, double del_y, double del_z)
{
	/*Hopping function*/
	double complex output = cexp(I * (q_x*del_x + q_y * del_y + q_z * del_z));
	return output;
}

double realdot(double *vec1, double *pos, int j)
{
	/*Function to take the dot product of two three-vectors*/
	double dot = 0;
	for(int i = 0; i<3; i++)
	{
		dot += vec1[i] * pos[3*j+i];
	}
	return dot;
}

int pos_idx(int i, int j)
{
	/* Function to index 1D arrays */
	int nn = 3*i + j;
	return nn;
}

void k_eff_field(double *EffFieldK, double *Dir, double *position, int N, double alpha, 
				double Dipol, int k_max)
{
	/* Function to calculate contribution to effective field from k-space */
	int k_x, k_y, k_z, m, n;	//Counting varibales
	for(m=0; m<(N*3); m++)		//Set array back to zero
	{
		EffFieldK[m] = 0;
	}
	/*Loop over allowed k vectors */
	for(k_x= -k_max; k_x <= k_max; k_x++)
	{
		for(k_y= -k_max; k_y <= k_max; k_y++)
		{
			for(k_z= -k_max; k_z <= k_max; k_z++)
			{
				/* Generate k vector and its magnitude */
				double kx = (double)k_x; double ky = (double)k_y; double kz = (double)k_z;
				double k_XYZ[3] = {kx, ky, kz};
				double k2 = realdot(k_XYZ, k_XYZ, 0);
				double k = sqrt(k2);
				/* Determine if the computed k-vector is allowed */
				if((k2<0.1) || (k2> (k_max*k_max +1)))
				{
					continue;
				}
				/* Compute k const outside loop */
				double kconst = 1.0/(8.0*8.0)*Dipol*pow(6.0, 1.5)*M_PI/k2*exp(-((M_PI*k)/(8.0*alpha))*((M_PI*k)/(8.0*alpha)));
				/* Loop over lattice. Parallelised for speed */
				#pragma omp parallel for default(shared) private(n) schedule(static, 1)
				for(m=0; m<N; m++)
				{
					for(n=0; n<N; n++)
					{
						double delx = position[3*m] - position[3*n];
						double dely = position[3*m+1] - position[3*n+1];
						double delz = position[3*m+2] - position[3*n+2];
						double Rij[3] = {delx, dely, delz};
						/*Add contributions to effective field */
						EffFieldK[3*m] += kconst * k_XYZ[0] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
						EffFieldK[3*m+1] += kconst * k_XYZ[1] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
						EffFieldK[3*m+2] += kconst * k_XYZ[2] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
					}
				}
			}
		}
	}
	return;
}

void local_k_eff_field(double *EffFieldK, double *Dir, double *position, int N, double alpha, 
				double Dipol, int k_max, int ion)
{
	/* Function to calculate contribution to effective field from k-space */
	int k_x, k_y, k_z, m, n;	//Counting varibales
	double *local_k_field = &EffFieldK[3*ion];
	for(m=0; m<3; m++)		//Set array back to zero
	{
		local_k_field[m] = 0;
	}
	/*Loop over allowed k vectors */
	for(k_x= -k_max; k_x <= k_max; k_x++)
	{
		for(k_y= -k_max; k_y <= k_max; k_y++)
		{
			for(k_z= -k_max; k_z <= k_max; k_z++)
			{
				/* Generate k vector and its magnitude */
				double kx = (double)k_x; double ky = (double)k_y; double kz = (double)k_z;
				double k_XYZ[3] = {kx, ky, kz};
				double k2 = realdot(k_XYZ, k_XYZ, 0);
				double k = sqrt(k2);
				/* Determine if the computed k-vector is allowed */
				if((k2<0.1) || (k2> (k_max*k_max +1)))
				{
					continue;
				}
				/* Compute k const outside loop */
				double kconst = 1.0/(8.0*8.0)*Dipol*pow(6.0, 1.5)*M_PI/k2*exp(-((M_PI*k)/(8.0*alpha))*((M_PI*k)/(8.0*alpha)));
				/* Loop over lattice. Parallelised for speed */
				for(n=0; n<N; n++)
				{
					double delx = position[3*ion] - position[3*n];
					double dely = position[3*ion+1] - position[3*n+1];
					double delz = position[3*ion+2] - position[3*n+2];
					double Rij[3] = {delx, dely, delz};
					/*Add contributions to effective field */
					local_k_field[0] += kconst * k_XYZ[0] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
					local_k_field[1] += kconst * k_XYZ[1] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
					local_k_field[2] += kconst * k_XYZ[2] * realdot(k_XYZ, Dir, n) * creal(cexp(M_PI*I*realdot(k_XYZ, Rij, 0)/4.0));
				}
			}
		}
	}
/*
	if(ion == 0)
	{
		printf("k space contribution to field at ion 0 is    (%lf, %lf, %lf)\n", EffFieldK[3*ion], EffFieldK[3*ion+1], EffFieldK[3*ion+2]);
	}
*/
	return;
}

void local_R_eff_field(double *EffFieldR, double *Dir, double *position, int N, double alpha, 
				double Dipol, int L_max, double J1, int NCells, int ion)
{
	const double J2 = -0.003;
	const double J3 = 0.013;
	int n_x, n_y, n_z, m, n;
	double *local_r_field = &EffFieldR[3*ion];
	for(m=0; m<3; m++)
	{
		local_r_field[m] = 0;
	}
	double six32 = sqrt(6*6*6);
	double cutoff = 64*L_max*L_max +1;
	double rconst = 2*Dipol*six32;
	for(n_x=-L_max; n_x<=L_max; n_x++)
	{
		for(n_y=-L_max; n_y<=L_max; n_y++)
		{
			for(n_z=-L_max; n_z<=L_max; n_z++)
			{
				double nx = (double)n_x; double ny = (double)n_y; double nz = (double)n_z;
				double delxyz[3] = {NCells*8*nx, NCells*8*ny, NCells*8*nz};
				for(n=0; n<N; n++)
				{
					double delx = position[3*ion] - position[3*n];
					double dely = position[3*ion+1] - position[3*n+1];
					double delz = position[3*ion+2] - position[3*n+2];
					double Rij[3] = {delx + delxyz[0], dely + delxyz[1], delz + delxyz[2]};
					double dist2 = realdot(Rij, Rij, 0);
					double distij = sqrt(dist2);
					if(5.99<dist2 && dist2<6.01)
					{
						EffFieldR[3*ion] += 2*J1*Dir[3*n];
						EffFieldR[3*ion+1] += 2*J1*Dir[3*n+1];
						EffFieldR[3*ion+2] += 2*J1*Dir[3*n+2];
					}
					if(13.99<dist2 && dist2<14.01)
					{	
						EffFieldR[3*ion] += 2*J2*Dir[3*n];
						EffFieldR[3*ion+1] += 2*J2*Dir[3*n+1];
						EffFieldR[3*ion+2] += 2*J2*Dir[3*n+2];
					}
					if(15.99<dist2 && dist2<16.01)
					{	
						EffFieldR[3*ion] += 2*J3*Dir[3*n];
						EffFieldR[3*ion+1] += 2*J3*Dir[3*n+1];
						EffFieldR[3*ion+2] += 2*J3*Dir[3*n+2];
					}
					if((0.1 < dist2) && (dist2 < cutoff))
					{
						// Add contribution to effective field
						EffFieldR[3*ion] += rconst*(BFun(distij, alpha)*Dir[3*n] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[0]);
						EffFieldR[3*ion+1] += rconst*(BFun(distij, alpha)*Dir[3*n+1] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[1]);
						EffFieldR[3*ion+2] += rconst*(BFun(distij, alpha)*Dir[3*n+2] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[2]);
					}
				}
			}
		}
	}
/*
	if(ion == 0)
	{
		printf("Real space contribution to field at ion 0 is (%lf, %lf, %lf)\n", EffFieldR[3*ion], EffFieldR[3*ion+1], EffFieldR[3*ion+2]);
	}
*/
	return;
}

void R_eff_field(double *EffFieldR, double *Dir, double *position, int N, double alpha, 
				double Dipol, int L_max, double J1, int NCells)
{
	int n_x, n_y, n_z, m, n;
	for(m=0; m<(N*3); m++)
	{
		EffFieldR[m] = 0;
	}
	double six32 = pow(6, 1.5);
	double cutoff = 64*L_max*L_max +1;
	double rconst = 2*Dipol*six32;
	for(n_x=-L_max; n_x<=L_max; n_x++)
	{
		for(n_y=-L_max; n_y<=L_max; n_y++)
		{
			for(n_z=-L_max; n_z<=L_max; n_z++)
			{
				double nx = (double)n_x; double ny = (double)n_y; double nz = (double)n_z;
				double delxyz[3] = {NCells*8*nx, NCells*8*ny, NCells*8*nz};
				#pragma omp parallel for default(shared) private(n) schedule(static, 1)
				for(m=0; m<N; m++)
				{
					for(n=0; n<N; n++)
					{
						double delx = position[3*m] - position[3*n];
						double dely = position[3*m+1] - position[3*n+1];
						double delz = position[3*m+2] - position[3*n+2];
						double Rij[3] = {delx + delxyz[0], dely + delxyz[1], delz + delxyz[2]};
						double dist2 = realdot(Rij, Rij, 0);
						double distij = sqrt(dist2);
						if(5.99<dist2 && dist2<6.01)
						{
							EffFieldR[3*m] += 2*J1*Dir[3*n];
							EffFieldR[3*m+1] += 2*J1*Dir[3*n+1];
							EffFieldR[3*m+2] += 2*J1*Dir[3*n+2];
						}
						if((0.1 < dist2) && (dist2 < cutoff))
						{
							// Add contribution to effective field
							EffFieldR[3*m] += rconst*(BFun(distij, alpha)*Dir[3*n] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[0]);
							EffFieldR[3*m+1] += rconst*(BFun(distij, alpha)*Dir[3*n+1] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[1]);
							EffFieldR[3*m+2] += rconst*(BFun(distij, alpha)*Dir[3*n+2] - CFun(distij, alpha)*realdot(Rij, Dir, n)*Rij[2]);
						}
					}
				}
			}
		}
	}
	return;
}

void Eff_field(double BB[3], double *EffField, double *EffFieldK, double *EffFieldR, 
				double *Dir, int N, double alpha, double Dipol, double Mag, double Zeeman)
{
	int i;
	for(i=0; i<(N*3); i++)
	{
		EffField[i] = 0;
	}
	double six32 = pow(6, 1.5);
	double FieldConst = Dipol*six32*(4*alpha*alpha*alpha)/(3*sqrt(M_PI));
	#pragma omp parallel for default(shared) schedule(static, 1)
	for(i=0; i<N; i++)
	{
		EffField[3*i] = -2*Mag*Zeeman*BB[0] + EffFieldK[3*i] + EffFieldR[3*i] - FieldConst*Dir[3*i];
		EffField[3*i+1] = -2*Mag*Zeeman*BB[1] + EffFieldK[3*i+1] + EffFieldR[3*i+1] - FieldConst*Dir[3*i+1];
		EffField[3*i+2] = -2*Mag*Zeeman*BB[2] + EffFieldK[3*i+2] + EffFieldR[3*i+2] - FieldConst*Dir[3*i+2];
	}
	return;
}

void local_Eff_field(double BB[3], double *EffField, double *EffFieldK, double *EffFieldR, 
				double *Dir, int N, double alpha, double Dipol, double Mag, double Zeeman, int ion)
{
	int i;
	double *local_eff_field = &EffField[3*ion];
	for(i=0; i<3; i++)
	{
		local_eff_field[i] = 0;
	}
	double six32 = sqrt(6*6*6);
	double FieldConst = Dipol*six32*(8*alpha*alpha*alpha)/(3*sqrt(M_PI));
	for(i=0; i<3; i++)
	{
		local_eff_field[0] = -2*Mag*Zeeman*BB[0] + EffFieldK[3*ion] + EffFieldR[3*ion] - FieldConst*Dir[3*ion];
		local_eff_field[1] = -2*Mag*Zeeman*BB[1] + EffFieldK[3*ion+1] + EffFieldR[3*ion+1] - FieldConst*Dir[3*ion+1];
		local_eff_field[2] = -2*Mag*Zeeman*BB[2] + EffFieldK[3*ion+2] + EffFieldR[3*ion+2] - FieldConst*Dir[3*ion+2];
	}
/*
	if(ion == 0)
	{
		printf("Total contribution to field at ion 0 is      (%lf, %lf, %lf)\n", EffField[3*ion], EffField[3*ion+1], EffField[3*ion+2]);
	}
*/
	return;
}

double Energy(double BB[3], double *EffField, double *Dir, double *energies, int N, double Mag, double Zeeman)
{
	int i, j;
	double E_tot=0;
	#pragma omp parallel for default(shared) private(j) reduction(+:E_tot) schedule(static, 1)
	for(i=0; i<N; i++)
	{
		double e_ion = 0;
		for(j=0; j<3; j++)
		{
			e_ion += EffField[3*i+j]*Dir[3*i+j] +2*Mag*Zeeman*BB[j]*Dir[3*i+j];
		}
		energies[i] = 0.5*e_ion;
		E_tot += energies[i];
	}
	return E_tot/N;
}

void normalise(double *field, int ion)
{
	double norm_sq = field[3*ion]*field[3*ion] + field[3*ion+1]*field[3*ion+1] + field[3*ion+2]*field[3*ion+2];
	double inv_norm = 1/sqrt(norm_sq);
	field[  3*ion  ] = inv_norm*field[  3*ion  ];
	field[3*ion + 1] = inv_norm*field[3*ion + 1];
	field[3*ion + 2] = inv_norm*field[3*ion + 2];
	return;
}

void cart_to_polar(double *cart, double polar[])
{
	double theta, phi, x, y, z, sign;
	x = cart[0];
	y = cart[1];
	z = cart[2];
	sign = x/fabs(x);
	theta = acos(z);
	phi = sign*atan(y/x);
	polar[0] = theta;
	polar[1] = phi;
	return;	
}

void polar_to_cart(double *cart, double polar[])
{
	double sint = sin(polar[0]);
	double cost = cos(polar[0]);
	double sinp = sin(polar[1]);
	double cosp = cos(polar[1]);
	
	cart[0] = cosp*sint;
	cart[1] = sint*sinp;
	cart[2] = cost;
	return;	
}

void Noise(double zhat[], double size)
{
	double d_sinth = size*genrand();
	double d_phi = 2*M_PI*genrand();
	double d_costh = sqrt(1- d_sinth*d_sinth);
	zhat[0] = cos(d_phi) * d_sinth;
	zhat[1] = sin(d_phi) * d_sinth;
	zhat[2] = d_costh;
}

void Perturb(double *Dir, int N, double size)
{
	int i;	//Counting variable
	for(i=0; i<N; i++)
	{
		double *spin = &Dir[3*i];
		double polar[2];
		double noise[3];	
	
		cart_to_polar(spin, polar);

		Noise(noise, size);

		//printf("noise is (%lf, %lf, %lf)\n", noise[0], noise[1], noise[2]);

		double sint = sin(polar[0]);
		double cost = cos(polar[0]);
		double sinp = sin(polar[1]);
		double cosp = cos(polar[1]);

		spin[0] = cost*cosp*noise[0] - sinp*noise[1] + sint*cosp*noise[2];
		spin[1] = cost*sinp*noise[0] - cosp*noise[1] + sint*sinp*noise[2];
		spin[2] = cost*noise[2] - sint*noise[0];

	}
	// Look to implement pertubation to spins here
	return;
}

void Cross(double *vec1, double vec2[], double out[])
{
	out[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1]; 
	out[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
	out[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
	return;
}

double dot(double *vec1, double *vec2)
{
	double out = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
	return out;
}
