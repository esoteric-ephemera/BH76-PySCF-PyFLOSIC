import yaml
from os import listdir,path,system,sys

"""
setup.py usage:
    python3 setup.py dfa=dfa_str clean=t/f

    dfa_str should either be a keyword in dict libxc_key (see below), or a
        LibXC format string, like GGA_X_PW91, MGGA_C_RSCAN
        This is the density functional approximation you want to use
        Defaults to None if not specified

    clean = boolean, only first character matters (t=True, f=False)
        Defaults to True

    examples : python3 setup.py dfa=LSDA clean=t
            python3 setup.py dfa=GGA_X_PW91,GGA_C_PW91
"""

base_dir = './'

libxc_key = {
    'LSDA': 'LDA_X, LDA_C_PW', 'PBE': 'GGA_X_PBE, GGA_C_PBE',
    'r2SCAN': 'MGGA_X_R2SCAN, MGGA_C_R2SCAN',
    'SCAN': 'MGGA_X_SCAN, MGGA_C_SCAN',
    'BLYP': 'GGA_X_B88, GGA_C_LYP'
}

# add eigenvalue level shifts here if you need them to stabilize convergence
lshiftd = { 'cl': 0.05 }

z_d = dict( H = 1, He = 2,
Li = 3, Be = 4, B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10,
Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18,
K = 19, Ca = 20, Sc = 21, Ti  = 22, V = 23, Cr = 24, Mn = 25, Fe = 26, Co = 27,
    Ni = 28, Cu = 29, Zn = 30, Ga = 31, Ge = 32, As = 33, Se = 34, Br = 25, Kr = 36,
Rb = 37, Sr = 38, Y = 39, Zr = 40, Nb = 41, Mo = 42, Tc = 43, Ru = 44, Rh = 45,
    Pd = 46, Ag = 47, Cd = 48, In = 49, Sn = 50, Sb = 51, Te = 52, I = 53, Xe = 54,
Pb = 82, Bi = 83 # I can't add all of them...
)

def get_ats_from_xyz(flnm):

    atd = {}
    with open(flnm,'r') as tfl:

        for iln,aln in enumerate(tfl):
            tmp = [x.strip() for x in aln.split()]
            if iln == 0:
                nion = int(aln.strip())
            elif 1 < iln < nion + 2:
                tmp = aln.strip().split()
                elt = tmp[0].strip()
                if len(elt) == 1:
                    elt = elt.upper()
                elif len(elt) == 2:
                    elt = elt[0].upper()+elt[1].lower()
                else:
                    print('Something is up, element name {:} too long or short'.format(elt))
                    raise SystemExit()
                if elt in atd:
                    atd[elt] += 1
                else:
                    atd[elt] = 1
    return atd


def setup(dfa,startclean=True):

    cmdict = yaml.load(open(base_dir+'BH76_chg_2s.yaml','r'), Loader=yaml.Loader)

    wdir = base_dir + '{:}_BH76_FODs/'.format(dfa)

    if startclean and path.isdir(wdir):
        system('rm -rf '+ wdir)

    """
    Specific options used universally:
        gridsize = 1,2,...,9 : grid density, 9 is most dense, 1 is least dense

        symmetry = boolean : use symmetry group of molecule to simplify calculation

        tol = float : convergence tolerance in Ha

        verbose = 0,1,...,9 : level of output from PySCF. 4 is a reasonable amout

        opt = FIRE, lbfgs, bfgs, linesearch, cg, gpmin : method for optimizing FODs

        max_inc = float : maximum change in FOD coordinates per iteration

        label = str : name of output files

        fmax = float : all |forces| < fmax for convergence

        max_step = integer : maximum number of steps ASE optimzation takes

        max_cycle = integer : maximum number of steps PySCF DFT calc takes
    """
    univ_opts = {
        'gridsize': 3, 'symm': False, 'tol': 1.e-7,
        'verbose': 4, 'opt': 'FIRE', 'max_inc' : 0.1, 'label': 'OPT_FRMORB',
        'fmax' : 1.e-4, 'max_step' : 1000, 'max_cycle': 500
    }

    srcdir = base_dir + 'BH76_geometries/'
    system('cp -r {:} {:}'.format(srcdir,wdir))
    system('cp sub_all_jobs.sh '+ wdir)

    for asys in listdir(srcdir):
        if asys[0] == '.':
            continue

        cdir = wdir + '/' + asys + '/'
        system('cp get_opt_fod.py '+ cdir)
        system('cp runjob.sh '+ cdir)

        optstr = ''
        for akey in univ_opts:
            optstr += '{:} = {:}\n'.format(akey,univ_opts[akey])

        optstr += 'charge = {:}\n'.format(cmdict[asys]['Charge'])
        optstr += '2S = {:}\n'.format(cmdict[asys]['2S'])

        if dfa in libxc_key:
            xcstr = libxc_key[dfa]
        else:
            xcstr = dfa
        optstr += 'xc = {:}\n'.format(xcstr)

        if asys in lshiftd:
            optstr += 'levelshift = {:}\n'.format(lshiftd[asys])

        optstr += 'xyz_fl = ./struc.xyz\n'

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
        basstr = 'DFO'
        optstr += 'basis = {:}\n'.format(basstr)

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

        #optstr += 'logfile = ' + asys + '.txt\n'

        inpfl = cdir + 'inp.txt'
        with open(inpfl,'w+') as tfile:
            tfile.write(optstr)

    return

def parse_into_setup():

    dfa = None
    makeclean = True
    if len(sys.argv) < 2:
        raise SystemExit('Need to specify density functional!')
    else:
        for astr in sys.argv[1:]:
            key,val = astr.split('=')
            if key.lower() == 'dfa':
                dfa = val
            elif key.lower() == 'clean':
                if val.lower()[0] == 't':
                    makeclean = True
                elif val.lower()[0] == 'f':
                    makeclean = False
                else:
                    print('Keyword "clean" must be t/f!')
            else:
                print('Skipping unknown keyword {:}'.format(key))

    setup(dfa,startclean=makeclean)

    return

if __name__ == "__main__":

    parse_into_setup()
