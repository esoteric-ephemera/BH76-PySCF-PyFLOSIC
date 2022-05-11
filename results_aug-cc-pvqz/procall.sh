#!/bin/bash
for u in ./*/ ; do
  cd $u
  python3 ../../analysis.py
  cd ..
done
