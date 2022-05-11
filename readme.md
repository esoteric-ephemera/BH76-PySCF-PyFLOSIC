## BH76 with PySCF and PyFLOSIC
#### Author: Aaron D. Kaplan (kaplan [@] temple.edu)

**Overview of files/folders**

+ BH76_chg_2s.yaml
  > YAML file containing all charge, multiplicity (2S, not 2S-1) and closed/open shell (R = restricted) data for the BH76 systems.
  Read in like a dictionary in Python (and human readable, unlike a pickle file)

  I had to rename the system c2h5 to c2h5_34 because on macOS, c2h5 and C2H5 are treated as equivalent directories (not so on Linux).

+ BH76_ref_energies.yaml
  > YAML file containing all reference energies and stoichiometries for the BH76 reactions

+ BH76_geometries
  > all .xyz and coord files for the BH76 set, taken from the GMTKN database:
  http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55/BH76.html

+ run_single_point.py
  > driver for running SCF DFT calculations with PySCF.
    + Input is parsed from inp.txt
    + Provided a calculation completes, a YAML file pyscf_run.yaml is written, containing the total energy and stating whether a calculation converged

+ run_single_point_hf.py
  > Same as run_single_point.py, but for HF calculation

+ setup.py
  > run this script to setup calculation for all systems in BH76 for a single density functional approximation.
  Examples are included in the file

+ analysis.py
  > run this in the current DFA directory to get all energies, errors, and error statistics

+ process.py
  > does a comparison between the DFT, DFT@HF, and DFT@HF(TS) (where HF density matrix is used only for transition states) errors

+ runjob.sh
  > sample job script to run all calculations within a directory, sequentially

+ results_aug-cc-pvqz
  > contains all data for aug-cc-pvqz basis set

+ FLOSIC
  > FOD optimization with PyFLOSIC. See readme in that directory for more info
