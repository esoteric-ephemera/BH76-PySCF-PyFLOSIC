#!/bin/bash
#PBS -l walltime=1:00:00:00
#PBS -q normal
#PBS -l nodes=1:ppn=20
#PBS -N BH76
#PBS -j oe
#PBS -o oe.txt
#PBS -m a

module --force purge ; module load intel-libs ; module load python/3.9.4
cd "$PBS_O_WORKDIR"

pyexec=/home/tuf53878/BH76/run_single_point.py

export OMP_NUM_THREADS=20
#export MKL_NUM_THREADS=20
#export OPENBLAS_NUM_THREADS=20
for u in ./*/ ; do
  cd $u
  python3 $pyexec
  cd ..
done
