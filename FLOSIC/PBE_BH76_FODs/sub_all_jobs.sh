#!/bin/bash

for u in ./*/ ; do
  cd $u
  qsub runjob.sh
  cd ..
done
