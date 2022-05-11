#!/bin/bash

srcs=(BH76_chg_2s.yaml get_opt_fod.py runjob.sh setup.py sub_all_jobs.sh BH76_geometries modified_pyflosic_files)

tar --disable-copyfile -cvf BH76_FOD_opt.tar ${srcs[*]}
