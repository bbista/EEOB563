#!/bin/bash


#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mail-user=bbista@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="raxmlrun.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="raxmlerrortest.txt" # job standard error file (%j replaced by job id)

module purge

module load intel

module load allinea

ulimit -s unlimited



module load raxml

raxmlHPC-PTHREADS-SSE3 -T2  -m GTRGAMMA -p 12345 -b 12345 -# 100 -s cob.aln -n T3
raxmlHPC-PTHREADS-SSE3 -T2  -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T1 -z RAxML_bootstrap.T3 -n T4
raxmlHPC-PTHREADS-SSE3 -T2  -m PROTGAMMAWAG -p 12345 -b 12345 -# 100 -s cobaa.aln -n T5
raxmlHPC-PTHREADS-SSE3 -T2 -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T2 -z RAxML_bootstrap.T5 -n T6