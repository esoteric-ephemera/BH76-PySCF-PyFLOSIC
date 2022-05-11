## Overview of files/folders:
#### Author: Aaron D. Kaplan (kaplan [@] temple.edu)

+ BH76_chg_2s.yaml
  > YAML file containing all charge, multiplicity (2S, not 2S-1) and closed/open shell (R = restricted) data for the BH76 systems.
  Read in like a dictionary in Python (and human readable, unlike a pickle file)

  I had to rename the system c2h5 to c2h5_34 because on macOS, c2h5 and C2H5 are treated as equivalent directories (not so on Linux).

+ BH76_geometries
  > all .xyz and coord files for the BH76 set, taken from the GMTKN database:
  http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/BH76.html

+ modified_pyflosic_files
  > A few files in PyFLOSIC
  (https://github.com/pyflosic/pyflosic)
  had to be modified:

  + preopt.py
    >version check fails if a version of PySCF > 1 is used.
    I fixed this so that any version of PySCF >= 1.5.2 passes

  + pycom.py
    > I've noticed that for some reason, the LSDA in PySCF needs a levelshift to converge for Cl.
    Added in the option to use a levelshift in the calculation of COM FODs.

+ get_opt_fod.py
  > this is the main driver that reads in a molecular geometry and calculates optimized FODs with PyFLOSIC.
  Computational parameters are set in inp.txt

  > If a run is successful, the center of mass (COM) FODs are written to \
  **FRMORB_COM** and **FB_GUESS_COM.xyz**,\
  > and the optimized FODs are written to \
  **FRMORB_OPT** and **FOD_OPT.xyz.xyz**

+ setup.py
  > run this script to setup FOD optimization for all systems in BH76 for a single density functional approximation.
  Examples are included in the file

+ runjob.sh
  > sample job script to submit. Modify this for your own use

+ sub_all_jobs.sh
  > simple file to submit all jobs in one DFA directory

+ LiH_LSDA_test
  > Example calculation with all input and output for LiH at lower convergence settings
