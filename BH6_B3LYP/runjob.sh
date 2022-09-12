#!/bin/bash
#PBS -l walltime=1:00:00:00
#PBS -q normal
#PBS -l nodes=1:ppn=20
#PBS -N BH6
#PBS -j oe
#PBS -o oe.txt
#PBS -m a

cd "$PBS_O_WORKDIR"
export PATH=/home/tuf53878/orca/openmpi-4.1.1/bin:$PATH
export LD_LIBRARY_PATH=/home/tuf53878/orca/openmpi-4.1.1/lib:$LD_LIBRARY_PATH
orca_exec=/home/tuf53878/orca/orca

for u in ./*/ ; do
  cd $u
  $orca_exec inp.txt >> out.txt
  cd ..
done