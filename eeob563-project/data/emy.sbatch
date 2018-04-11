#!/bin/bash


#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mail-user=bbista@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="codeout.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="codeerror.txt" # job standard error file (%j replaced by job id)

module purge

module load intel

module load allinea

ulimit -s unlimited



module load R/3.3.1 

module load prank/150803

module load paml/4.9

Rscript code.R

