#!/bin/bash

srcs=(BH76_chg_2s.yaml BH76_ref_energies.yaml run_single_point.py runjob.sh setup.py analysis.py BH76_geometries)

tar --disable-copyfile -cvf bh76_pyscf.tar ${srcs[*]}
