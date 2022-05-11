import yaml

from os import listdir,path,system,sys

base_dir = './'
hf_wvfn_dir = './HF_WVFN_aug-cc-pvqz'

xclib = 'LibXC'

xcfun_key = {
    'LSDA': '1.0*SLATERX, 1.0*PW92C', 'PBE': '1.0*PBEX, 1.0*PBEC',
    'r2SCAN': '1.0*R2SCANX, 1.0*R2SCANC', 'SCAN': '1.0*SCANX, 1.0*SCANC',
    'BLYP': '1.0*BECKEX, 1.0*LYPC'
}

libxc_key = {
    'LSDA': 'LDA_X, LDA_C_PW', 'PBE': 'GGA_X_PBE, GGA_C_PBE',
    'r2SCAN': 'MGGA_X_R2SCAN, MGGA_C_R2SCAN',
    'SCAN': 'MGGA_X_SCAN, MGGA_C_SCAN',
    'BLYP': 'GGA_X_B88, GGA_C_LYP'
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

    wdir = base_dir + '{:}_BH76/'.format(dfa)

    if startclean and path.isdir(wdir):
        system('rm -rf '+ wdir)

    # specific options used universally
    univ_opts = {
        'gridsize': 9, 'symm': False, 'tol': 1.e-7, 'max_cycle': 500,
        'verbose': 4
    }

    srcdir = base_dir + 'BH76_geometries/'
    system('cp -r {:} {:}'.format(srcdir,wdir))
    system('cp runjob.sh '+ wdir)
    if dfa == 'HF':
        system('mkdir {:}HF_WVFN'.format(wdir))

    hfdft = False
    if dfa[-3:] == '@HF':
        dfa = dfa[:-3]
        univ_opts['max_cycle'] = 0
        hfdft = True
        if path.isdir(hf_wvfn_dir) and len(listdir(hf_wvfn_dir)) > 0:
            for asys in cmdict:
                if asys == 'h':
                    # we use reduced density matrix for h, PySCF won't store the
                    # wavefunction file here. Dunno why
                    continue
                cstr = 'cp ./{:}/{:}_wvfn {:}{:}/wvfn'.format(hf_wvfn_dir,asys,\
                    wdir,asys)
                system(cstr)
        else:
            raise SystemExit("I can't do DFA@HF without Hartree-Fock orbitals!")

    for asys in listdir(srcdir):
        if asys[0] == '.':
            continue

        cdir = wdir + '/' + asys + '/'
        if dfa == 'HF':
            system('cp run_single_point_hf.py '+ cdir + 'run_single_point.py')
            univ_opts.pop('gridsize',None)
        else:
            system('cp run_single_point.py '+ cdir)

        optstr = ''
        for akey in univ_opts:
            optstr += '{:} = {:}\n'.format(akey,univ_opts[akey])

        optstr += 'charge = {:}\n'.format(cmdict[asys]['Charge'])
        optstr += '2S = {:}\n'.format(cmdict[asys]['2S'])
        if cmdict[asys]['R']:
            optstr += 'restricted = True\n'
        else:
            optstr += 'restricted = False\n'

        if dfa != 'HF':
            if xclib == 'XCFun':
                optstr += 'xc = {:}\n'.format(xcfun_key[dfa])
                optstr += 'xc_lib = XCFun\n'
            elif xclib == 'LibXC':
                optstr += 'xc = {:}\n'.format(libxc_key[dfa])
                optstr += 'xc_lib = LibXC\n'
            else:
                raise SystemExit('Unknown XC library '+xclib)

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
        basstr = 'aug-cc-pvqz'#'def2-QZVP'
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

        optstr += 'logfile = ' + asys + '.txt\n'
        if dfa == 'HF':
            optstr += 'checkfile = ../HF_WVFN/{:}_wvfn\n'.format(asys)

        if hfdft:
            if asys == 'h':
                optstr += 'init = HF DM\n'
            else:
                optstr += 'init = chkfl\n'
                optstr += 'chkfl = ./wvfn'

        inpfl = cdir + 'inp.txt'
        with open(inpfl,'w+') as tfile:
            tfile.write(optstr)

    return

if __name__ == "__main__":

    setup('LSDA@HF',startclean=True)
