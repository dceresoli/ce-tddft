#!/bin/bash --login
#======================================================================
#SBATCH --job-name=WAVEPACKET
#SBATCH --nodes=1
#SBATCH --ntasks=8 
#SBATCH --time=01:30:00
#SBATCH --partition=parallel,samsung,calcolo
#======================================================================
scratch=/scratch/${USER}${PWD#$HOME}
mkdir -p $scratch 2>/dev/null
ln -sf $scratch ./scratch
#======================================================================
module purge 
module load intel/xe2015 openmpi/1.8.3-intel
export PATH=/home/ceresoli/Codes/q-e/bin:$PATH
export PATH=/home/ceresoli/Codes/ce-tddft/bin:$PATH
#======================================================================
mpirun pw.x -ndiag 1 <benzene-scf.in >benzene-scf.out
mpirun tddft.x -ndiag 1 <benzene-wp.in >benzene-wp.out


