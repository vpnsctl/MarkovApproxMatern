#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=125
#SBATCH --time=60:00:00


export OMP_NUM_THREADS=125

 module load gcc/12.2.0 
 module load fftw/3.3.10/openmpi-4.1.4-gcc-12.2.0-dp
 module load udunits/2.2.28 
 module load openblas/0.3.21/gnu-12.2.0 
 module load curl/7.86.0/gnu-12.2.0
 module load R/4.3.2/gnu-12.2.0
 module load openmpi/4.1.4/gnu11.2.1 
 module load binutils/2.37
 module load geos/3.12.1 
 module load gsl/2.7.1/gnu-12.2.0 
 module load gdal/3.6.2/gnu-12.2.0  
 
srun -c 125 R CMD BATCH --vanilla compute_pred_error_5000_fourier_range1.R