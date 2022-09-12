#!/bin/bash
#PBS -l walltime=1:00:00:00
#PBS -q normal
#PBS -l nodes=1:ppn=20
#PBS -N BH76_FOD_opt
#PBS -j oe
#PBS -o oe.txt
#PBS -m a

module --force purge ; module load intel-libs ; module load python/3.9.4
cd "$PBS_O_WORKDIR"

export OMP_NUM_THREADS=20
#export MKL_NUM_THREADS=20
#export OPENBLAS_NUM_THREADS=20

pyfd=/home/tuf53878/pyflosic

export PYTHONPATH=$pyfd/src/:$PYTHONPATH
export PYTHONPATH=$pyfd/utils/:$PYTHONPATH

python3 run_sp.py
