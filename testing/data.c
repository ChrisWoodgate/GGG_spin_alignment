/*==============================================/
/                                               /						
/ Data handling routines for GGG spin alignment /
/ code.                                         /
/ Original code by C Woodgate - 2018            /
/												/
/==============================================*/

#include "data.h"

void read_input(int *N_ions, double *alpha, int *k_max, int *L_max, double *theta, double *phi)  
{
	//Function to read the input file.
	FILE *infile;
	if(!(infile=fopen("input.txt","r"))) 
	{   
		printf("Error opening file\n");
		exit(1);
	}   
	if(6!=fscanf(infile,"%d %lf %d %d %lf %lf", 
				N_ions, alpha, k_max, L_max, theta, phi)) 
	{   
		printf("Error reading parameters from file\n");
		exit(1);
	}   
	fclose(infile);
	return;
}

void read_positions(long N_ions, double *positions, int NCells) 
{
        //Function to read positions file.
        FILE *infile;
    	char filename[sizeof "positions/positions_00Cells.txt"];
    	sprintf(filename,"positions/positions_%02dCells.txt",NCells);
		if(!(infile=fopen(filename,"r"))) 
        {   
                printf("Error opening file\n");
                exit(1);
        }   
        for(int i=0; i<3*N_ions; i++)
        {   
                if(1!=fscanf(infile, "%lf \n", &(positions[i])))
                {   
                        printf("Error reading coefficients from file\n");
                        exit(1);
                }   
        }   
        fclose(infile);
		return;
}

/* Parse Command line arguments */
int read_args(int argc, char **argv, int *spin_movement, int *pertubation, int *NSteps, int *NCells, double *field_strength)
{
  int index;
  int c;
  opterr = 0;

  /* Process all flags found in argc */
  while ((c = getopt(argc, argv, "B:C:N:PI")) != -1)
    switch (c)
      {
      case 'I':
        *spin_movement = 1;
        break;
      case 'P':
        *pertubation = 1;
        break;
      case 'N':
        *NSteps = atoi(optarg);
        if (*NSteps<=0) {
           opterr = 1; 
           fprintf (stderr, "Number of steps argument could not be read: %s\n", optarg);
        }
        break;
      case 'C':
        *NCells = atoi(optarg);
        if (*NCells<=0) {
           opterr = 1; 
           fprintf (stderr, "Number of cells argument could not be read: %s\n", optarg);
        }
        break;
      case 'B':
        *field_strength = atof(optarg);
        if (*field_strength<0.0) {
           opterr = 1; 
           fprintf (stderr, "Field Strength argument could not be read: %s\n", optarg);
        }
        break;
      case '?':
        if ((optopt=='N')||(optopt=='B')||(optopt=='C'))
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
        return 1;
      default:
        abort ();
      }

  /* List unrecognised arguments */
  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);

  return opterr;
}
