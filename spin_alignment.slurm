#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2012mb
#SBATCH --time=24:00:00

# Change into the directory where I submitted the job
cd $SLURM_SUBMIT_DIR 

export OMP_NUM_THREADS=16

for p in `seq 3000 -100 1600`
do
	time ./GGG_spins -N 10 -C 5 -B $p
done
exit 0
