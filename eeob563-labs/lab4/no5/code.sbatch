#!/bin/bash


#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=32   # number of nodes
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

rm *.T3
rm *.T4
rm *.T5
rm *.T6

mpirun -np 16 /shared/class/eeob563/raxml/standard-RAxML/raxmlHPC-MPI  -m GTRGAMMA -p 12345 -b 12345 -# 100 -s cob.aln -n TT1
mpirun -np 16 /shared/class/eeob563/raxml/standard-RAxML/raxmlHPC-MPI  -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T1 -z RAxML_bootstrap.TT1 -n TT2
mpirun -np 16 /shared/class/eeob563/raxml/standard-RAxML/raxmlHPC-MPI  -m PROTGAMMAWAG -p 12345 -b 12345 -# 100 -s cobaa.aln -n TT5
mpirun -np 16 /shared/class/eeob563/raxml/standard-RAxML/raxmlHPC-MPI  -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T2 -z RAxML_bootstrap.TT5 -n TT6