from ase.units import Bohr
import yaml
from os import path, system, sys

z_d = dict( H = 1, He = 2,
    Li = 3, Be = 4, B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10,
    Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18,
    K = 19, Ca = 20, Sc = 21, Ti  = 22, V = 23, Cr = 24, Mn = 25, Fe = 26, Co = 27,
        Ni = 28, Cu = 29, Zn = 30, Ga = 31, Ge = 32, As = 33, Se = 34, Br = 25, Kr = 36,
    Rb = 37, Sr = 38, Y = 39, Zr = 40, Nb = 41, Mo = 42, Tc = 43, Ru = 44, Rh = 45,
        Pd = 46, Ag = 47, Cd = 48, In = 49, Sn = 50, Sb = 51, Te = 52, I = 53, Xe = 54,
    Pb = 82, Bi = 83 # I can't add all of them...
)

# just the directory where this file and BH76_chg_2s.yaml are located
base_dir = path.dirname(path.realpath(__file__)) + '/'

def read_xyz_file(flnm):

    outstr = ''
    with open(flnm,'r') as infl:
        for irow,arow in enumerate(infl):
            # irow = line # starting from 0
            # arow = contents of the line
            # strip() removes whitespace and newline character \n
            tval = arow.strip()
            if irow == 0:
                Natom = int(tval)
                # line 2 (irow=1) is comment line for xyz file
            elif irow > 1:
                # split(char) separates a string around character char
                tmpv = [x.strip() for x in tval.split()]
                elt = '' # element name
                pos = [float(x)/Bohr for x in tmpv[1:]] # positions in Bohr, not Angstrom!
                for ichar,achar in enumerate(tmpv[0]):
                    if ichar == 0:
                        elt += achar.upper()
                    else:
                        elt += achar.lower()
                for i in range(3):
                    if pos[i] < 0:
                        pad = ' '
                    else:
                        pad = '  '
                    outstr += '{:}{:.12f}'.format(pad,pos[i])

                zval = z_d[elt]
                if zval >= 10:
                    pad = ' '
                else:
                    pad = '  '
                outstr += '{:}{:} ALL\n'.format(pad,zval)
                #outstr += '{:.12f}  {:.12f}  {:.12f}  {:} ALL\n'.format(*pos,z_d[elt])
                #if elt not in at_l:
                #    at_l.append(elt)
    return outstr, Natom

def make_nrlmol_inp(uopts={},fname='./NRLMOL_INPUT.DAT'):

    defopts = dict(
        ATOMSPHV = 'N', BASISV = 'DEFAULT', CALCTYPEV = 'SCF-ONLY', DFTD3V = 'N',
        DIAG1V =  1, DIAG2V = 1, DIAG3V =  0, DMATV = 'N', DOSOCCUV = 'N',
        FIXMV = 'N', FOD_LOOPV = 'N', FOD_OPT1V = 'CG', FOD_OPT2V = 'N',
        FOD_OPT3V = 0, FRAGMENTV = 'N', JNTDOSV = 'N', MAXSCFV = 100, MIXINGV = 'P',
        NONSCFV = 'N', NONSCFFORCESV = 'N', NWFOUTV = 10, POPULATIONV = 'N',
        RHOGRIDV = 'N', SCFTOLV = 1.e-6, SPNPOLV = 'N', MESHSETTINGV = 1,
        SYMMETRYV = 'N', UNIHAMV = 'N', WFGRIDV = 'N', WFFRMV = 'N'
    )

    for akey in uopts:
        # replace or append default options with user specified options
        defopts[akey] = uopts[akey]

    mklen = 0
    for akey in defopts:
        mklen = max(mklen,len(akey))

    minpad = 2

    tstr = '#\n#\n#\n\n&input_data'
    for akey in defopts:
        tmpstr = '\n{:}'.format(akey) + ' '*(mklen-len(akey)+minpad)
        if isinstance(defopts[akey], str):
            tmpstr += "= '{:}'".format(defopts[akey])
        elif isinstance(defopts[akey], int):
            tmpstr += '= {:}'.format(defopts[akey])
        elif isinstance(defopts[akey], float):
            tmpstr += '= {:.1e}'.format(defopts[akey]).replace('e','D')

        tstr += tmpstr

    tstr += '\n&end\n\n'

    with open(fname,'w+') as tfl:
        tfl.write(tstr)

    return

def make_cluster(dfa_x,dfa_c,xyz_fl,chrg,spin,symm=None,fname='./CLUSTER'):
    """
    Function make_cluster
    Inputs:
        dfa_x (str) : exchange density functional
        dfa_c (str) : correlation density functional
        xyz_fl (str) : name of xyz file to be converted
        chrg (float) : charge of molecule
        spin (float) : 2*S_tot
        symm (str) : Optional, defaults to None. Symmetry of molecule

    Output:
        CLUSTER file
    """

    outstr = '{:}*{:}\n'.format(dfa_x,dfa_c)
    outstr += '{:}\n'.format(symm)

    geo_str, Natom = read_xyz_file(xyz_fl)
    outstr += '{:}\n'.format(Natom) # number of inequivalent atoms
    outstr += geo_str # geometry of the molecule
    outstr += '{:}  {:}\n'.format(chrg,spin)

    # write the CLUSTER file
    with open(fname,'w+') as tfl:
        tfl.write(outstr)

    return

def make_job_fl(dfa,sys,fname='./dft.sh'):
    job_fl_str = '#!/bin/bash\n'
    job_fl_str += '#PBS -l walltime=04:30:00\n'
    job_fl_str += '#PBS -q normal\n'
    job_fl_str += '#PBS -l nodes=1:ppn=20\n'
    job_fl_str += '#PBS -N {:}-{:}\n'.format(dfa,sys)
    job_fl_str += '#PBS -A cst\n\n'
    job_fl_str += 'flosic_exe=/home/tuj67482/flosic-aaron/FLOSIC_public_2020_r2SCAN/flosic/nrlmol_exe\n'
    job_fl_str += 'cd "$PBS_O_WORKDIR"\n'
    job_fl_str += 'mpirun $flosic_exe > output'

    with open(fname,'w+') as tfl:
        tfl.write(job_fl_str)

    return

def make_all_cluster_files(dfa,clean=True):

    dfad = {'LSDA': 'LDA-PW91', 'PBE': 'GGA-PBE', 'SCAN': 'MGGA-SCAN'}

    cmdict = yaml.load(open(base_dir+'BH76_chg_2s.yaml','r'), Loader=yaml.Loader)

    if dfa in dfad:
        # Use preformatted DFA, like LSDA
        dfax = dfad[dfa]
        dfac = dfad[dfa]
    else:
        # Use a DFA of the form DFA_X*DFA_C
        dfax = dfa.split('*')[0]
        dfac = dfa.split('*')[1]

    # {:} turns a variable into a string using format
    wdir = '{:}/BH76_{:}'.format(base_dir,dfa)
    if clean:
        system('rm -rf '+ wdir)

    for asys in cmdict:
        cdir = '{:}/{:}'.format(wdir,asys)
        system('mkdir -p '+cdir)
        chg = cmdict[asys]['Charge']
        two_s = cmdict[asys]['2S']
        filename = cdir+'/CLUSTER'
        geofile = base_dir+'/BH76_geometries/{:}/struc.xyz'.format(asys)
        # first make CLUSTER file
        make_cluster(dfax,dfac,geofile,chg,two_s,symm=None,fname=filename)
        # now make NRLMOL_INPUT.DAT file
        make_nrlmol_inp(fname=cdir+'/NRLMOL_INPUT.DAT')
        # and TMPTRE
        with open(cdir + '/TMPTRE','w+') as tfl:
            tfl.write('   9.9999999999999995E-007  KT IN HARTREES')
        # finally the job file
        make_job_fl(dfa,asys,fname=cdir+'/dft.sh')

    tstr = '#!/bin/bash\n# run me using\n#    bash sub_all_jobs.sh\n'
    tstr += 'for u in ./*/ ; do\n'
    tstr += '  cd $u \n  qsub dft.sh\n  cd ..\n'
    tstr += 'done'
    with open(wdir+'/sub_all_jobs.sh','w+') as tfl:
        tfl.write(tstr)

    return

if __name__ == "__main__":

    if len(sys.argv) < 2:
        raise SystemExit('WARNING: need to specify DFA')
    else:
        dfa = sys.argv[1]
    make_all_cluster_files(dfa)
    #make_cluster('GGA-PBE','LDA-PW91','./BH76_geometries/c2h4/struc.xyz',0.0,0.0)
