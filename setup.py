import yaml
from os import listdir,path,system,sys

from run_single_point import opt_type,str_parse,get_ats_from_xyz

base_dir = path.dirname(path.realpath(__file__)) + '/'

"""
setup.py usage:
    python3 setup.py dfa=str clean=t/f at_dfa=str dfa_dir=str

    dfa should either be a keyword in dict libxc_key (see below), or a
        LibXC format string, like GGA_X_PW91,MGGA_C_RSCAN
        This is the density functional approximation you want to use
        Defaults to None if not specified

    clean = boolean, only first character matters (t=True, f=False)
        Defaults to True

    xclib = LibXC or XCFun, sets library used to get XC DFA
        Defaults to LibXC

    at_dfa should be a directory name where wavefunction files are written to
        Used to do dfa@at_dfta calculations

    dfa_dir should be the name of the DFA directory if you don't want to use dfa as its name
        Defaults to be same as dfa

    examples : python3 setup.py dfa=LSDA clean=t
            python3 setup.py dfa=GGA_X_PW91,GGA_C_PW91 xclib=LibXC
            python3 setup.py dfa=PW91X,PW91C xclib=XCFun
"""

xcfun_key = {
    'LSDA': '1.0*SLATERX, 1.0*PW92C', 'PBE': '1.0*PBEX, 1.0*PBEC',
    'r2SCAN': '1.0*R2SCANX, 1.0*R2SCANC', 'SCAN': '1.0*SCANX, 1.0*SCANC',
    'BLYP': '1.0*BECKEX, 1.0*LYPC'
}

libxc_key = {
    'LSDA': 'LDA_X, LDA_C_PW', 'PBE': 'GGA_X_PBE, GGA_C_PBE',
    'r2SCAN': 'MGGA_X_R2SCAN, MGGA_C_R2SCAN',
    'SCAN': 'MGGA_X_SCAN, MGGA_C_SCAN',
    'BLYP': 'GGA_X_B88, GGA_C_LYP',
    'B3LYP': 'HYB_GGA_XC_B3LYP',
    'LCwPBE': 'HYB_GGA_XC_LC_WPBE',
    'M06-L': 'MGGA_X_M06_L, MGGA_C_M06_L',
    'MN15-L': 'MGGA_X_MN15_L, MGGA_C_MN15_L',
}

# add eigenvalue level shifts here if you need them to stabilize convergence
lshiftd = { }

z_d = dict( H = 1, He = 2,
Li = 3, Be = 4, B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10,
Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18,
K = 19, Ca = 20, Sc = 21, Ti  = 22, V = 23, Cr = 24, Mn = 25, Fe = 26, Co = 27,
Ni = 28, Cu = 29, Zn = 30, Ga = 31, Ge = 32, As = 33, Se = 34, Br = 25, Kr = 36,
Rb = 37, Sr = 38, Y = 39, Zr = 40, Nb = 41, Mo = 42, Tc = 43, Ru = 44, Rh = 45,
Pd = 46, Ag = 47, Cd = 48, In = 49, Sn = 50, Sb = 51, Te = 52, I = 53, Xe = 54,
Pb = 82, Bi = 83 # I can't add all of them...
)


def setup(dfa,startclean=True,xclib='LibXC',wfdir = '',dfa_dir=None,inp_opts={}):

    cmdict = yaml.load(open(base_dir+'BH76_chg_2s.yaml','r'), Loader=yaml.Loader)

    if dfa_dir is not None:
        wdir = base_dir + '{:}_BH76/'.format(dfa_dir)
    else:
        wdir = base_dir + '{:}_BH76/'.format(dfa)

    if startclean and path.isdir(wdir):
        system('rm -rf '+ wdir)

    # specific options used universally
    univ_opts = {
        'gridsize': 9, 'symm': False, 'tol': 1.e-7, 'max_cycle': 500,
        'verbose': 4, 'basis': 'aug-cc-pvqz'
    }
    for akey in inp_opts:
        univ_opts[akey] = inp_opts[akey]

    srcdir = base_dir + 'BH76_geometries/'
    system('cp -r {:} {:}'.format(srcdir,wdir))
    system('cp runjob.sh '+ wdir)
    if dfa == 'HF':
        system('mkdir {:}HF_WVFN'.format(wdir))

    nscf_dft = False
    if '@' in dfa:
        nscf_dft = True
        dfa,at_dfa = dfa.split('@')
        univ_opts['max_cycle'] = 0
        if path.isdir(wfdir) and len(listdir(wfdir)) > 0:
            for asys in cmdict:
                if asys == 'h' and at_dfa == 'HF':
                    # we use reduced density matrix for h, PySCF won't store the
                    # HF wavefunction file here. Dunno why
                    continue
                cstr = 'cp {:}{:}/{:}_wvfn {:}{:}/wvfn'.format(base_dir,wfdir,\
                    asys,wdir,asys)
                system(cstr)
        else:
            err_str = "I can't do {:}@{:} without {:} orbitals!".format(dfa,at_dfa,at_dfa)
            raise SystemExit(err_str)

    for asys in listdir(srcdir):
        if asys[0] == '.':
            continue

        cdir = wdir + '/' + asys + '/'
        #if dfa == 'HF':
        #    system('cp run_single_point_hf.py '+ cdir + 'run_single_point.py')
        #    univ_opts.pop('gridsize',None)
        #else:
        #    system('cp run_single_point.py '+ cdir)

        optstr = ''
        for akey in univ_opts:
            optstr += '{:} = {:}\n'.format(akey,univ_opts[akey])

        optstr += 'charge = {:}\n'.format(cmdict[asys]['Charge'])
        optstr += '2S = {:}\n'.format(cmdict[asys]['2S'])
        if cmdict[asys]['R']:
            optstr += 'restricted = True\n'
        else:
            optstr += 'restricted = False\n'

        if dfa == 'HF':
            optstr += 'xc = HF\n'
        else:
            xcstr = dfa
            if xclib == 'XCFun':
                if dfa in xcfun_key:
                    xcstr = xcfun_key[dfa]
            elif xclib == 'LibXC':
                if dfa in libxc_key:
                    xcstr = libxc_key[dfa]
            else:
                raise SystemExit('Unknown XC library '+xclib)

            optstr += 'xc = {:}\n'.format(xcstr)
            optstr += 'xc_lib = {:}\n'.format(xclib)

        if asys in lshiftd:
            optstr += 'levelshift = {:}\n'.format(lshiftd[asys])

        optstr += 'xyz_fl = ./struc.xyz\n'# + cdir + 'struc.xyz\n'

        strucfl = srcdir + '{:}/struc.xyz'.format(asys)
        chem_comp = get_ats_from_xyz(strucfl)
        #nats = len(chem_comp.keys())

        """
            you can modify the basis you want to use here
            if, for example, you want to use different bases for different atoms,
            use the format
                basis = atom 1 : basis 1 ; atom 2 : basis 2 ; ...
            Example: H2O. To use STO-3G with H and cc-pVDZ with O, write
                basis = H : STO-3G ; O : cc-pVDZ
            spaces are ignored, use enough to be readable
        """
        #basstr = 'aug-cc-pvqz'#'def2-QZVP'
        #optstr += 'basis = {:}\n'.format(basis)

        """
            For heavy atoms, an electronic core potential (ECP)/pseudopotential is
            used by default when the atomic Z > 36 (Kr).

            These aren't used in BH76 at all, but included here as an example
        """
        ecpstr = ''
        iat = 0
        mz = -1
        for anat in chem_comp:

            if iat == 0:
                cchar = ''
            else:
                cchar = '; '

            if z_d[anat] > z_d['Kr']:
                ecpstr += '{:}{:}: def2-qzvp'.format(cchar,anat)
                iat += 1

        if len(ecpstr) > 0:
            optstr += 'ecp = {:}\n'.format(ecpstr)

        optstr += 'logfile = ' + asys + '.txt\n'
        if dfa == 'HF':
            optstr += 'write_chkfl = True\n'
            optstr += 'chkfl = ../HF_WVFN/{:}_wvfn\n'.format(asys)

        if nscf_dft:
            if asys == 'h' and at_dfa == 'HF':
                optstr += 'init = HF DM\n'
            else:
                optstr += 'init = chkfl\n'
                optstr += 'chkfl = ./wvfn'

        inpfl = cdir + 'inp.txt'
        with open(inpfl,'w+') as tfile:
            tfile.write(optstr)

    return

def parse_into_setup():

    dfa = None
    makeclean = True
    xclib = 'LibXC'
    wfdir = ''
    dfa_dir = None

    ioptd = {}
    ioptl = []
    for akey in opt_type:
        for bkey in opt_type[akey]:
            ioptl.append(bkey)

    if len(sys.argv) < 2:
        raise SystemExit('Need to specify density functional!')
    else:
        for astr in sys.argv[1:]:
            key,val = astr.split('=')
            lkey = key.lower()
            if lkey == 'dfa':
                dfa = val
            elif lkey == 'clean':
                if val.lower()[0] == 't':
                    makeclean = True
                elif val.lower()[0] == 'f':
                    makeclean = False
                else:
                    print('Keyword "clean" must be t/f!')
            elif lkey == 'xclib':
                xclib = val
            elif lkey == 'at_dfa':
                wfdir = val
            elif lkey == 'dfa_dir':
                dfa_dir = val
            elif key in ioptl:
                ioptd[key] = val
            else:
                print('Skipping unknown keyword {:}'.format(key))

    setup(dfa,startclean=makeclean,xclib=xclib,wfdir=wfdir,dfa_dir=dfa_dir,\
        inp_opts=ioptd)

    return

if __name__ == "__main__":

    #setup('LSDA@HF',startclean=True)
    parse_into_setup()
