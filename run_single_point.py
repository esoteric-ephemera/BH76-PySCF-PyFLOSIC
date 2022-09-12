import numpy as np
from pyscf import gto,dft
from os import path

opt_type = {
    'float': ['tol','levelshift'],
    'int': ['gridsize','max_cycle','charge','2S','verbose'],
    'bool': ['symm','restricted','write_chkfl','cartesian'],
    'str': ['xyz_fl','basis','xc','xc_lib','ecp','logfile','init','chkfl']
}

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

def get_nrlmol_bas(xyz_file,var='dfo'):

    if var.lower() in ['nrlmol','dfo']:
        wfl = path.dirname(path.realpath(__file__))+'/dfo-nrlmol.dat'
    elif var.lower() in ['nrlmol+','dfo+']:
        wfl = path.dirname(path.realpath(__file__))+'/dfo+-nrlmol.dat'
    else:
        raise ValueError('Unknown basis set {:}; stopping.'.format(var))

    atd = get_ats_from_xyz(xyz_file)
    basd = {}
    for at in atd:
        basd[at] = gto.load(wfl,symb=at)
    return basd

def str_parse(key,val):
    if key in opt_type['str']:
        # quick return
        return val
    elif key in opt_type['float']:
        return float(val)
    elif key in opt_type['int']:
        return int(val)
    elif key in opt_type['bool']:
        if val.lower() in ['true','t']:
            return True
        elif val.lower() in ['false','f']:
            return False
        else:
            estr = "Could not process key {:} with value {:}; expected boolean".format(key,val)
            raise SystemExit(estr)
    else:
        estr = "Unknown key {:} with value {:}".format(key,val)
        raise SystemExit(estr)

def parse_inp(flnm):

    opts = {}
    with open(flnm,'r') as tfl:
        for row in tfl:
            tmp = [x.strip() for x in row.split('=')]
            if len(tmp) == 1:
                # blank line
                continue
            if ':' in tmp[1]:
                # dict type option
                tmp2 = [x.strip() for x in tmp[1].split(';')]
                opts[tmp[0]] = {}
                for x in tmp2:
                    tmp3 = [y.strip() for y in x.split(':')]
                    opts[tmp[0]][tmp3[0].strip()] = str_parse(tmp[0],tmp3[1])
            else:
                opts[tmp[0]] = str_parse(tmp[0],tmp[1])

    return opts

def molscf():

    dopts = {
        'gridsize': 5, 'basis': 'def2-QZVP', 'symm': False,
        'tol': 1.e-7, 'max_cycle': 500, 'charge': 0, '2S': 0,
        'xc': '1.0*SLATERX, 1.0*PW92C', 'xc_lib': 'XCFun', 'verbose': 4,
        'restricted': False, 'ecp' : {}, 'init': None, 'write_chkfl' : False,
        'cartesian': False
        }

    uopts = parse_inp('./inp.txt')

    for akey in uopts:
        dopts[akey] = uopts[akey]

    if dopts['basis'].lower() in ['nrlmol','dfo','nrlmol+','dfo+']:
        # NRLMOL basis sets
        dopts['basis'] = get_nrlmol_bas(dopts['xyz_fl'],var=dopts['basis'])

    mol = gto.M(atom=dopts['xyz_fl'], basis=dopts['basis'], symmetry=dopts['symm'], \
        charge=dopts['charge'], spin=dopts['2S'], output = dopts['logfile'], \
        verbose = dopts['verbose'], ecp = dopts['ecp'], cart = dopts['cartesian'])

    if dopts['xc'] == 'HF':
        from pyscf import scf
        if dopts['restricted']:
            kscalc = scf.RHF(mol)
        else:
            kscalc = scf.UHF(mol)
    else:
        if dopts['restricted']:
            kscalc = dft.RKS(mol)
        else:
            kscalc = dft.UKS(mol)

    if dopts['write_chkfl']:
        if 'chkfl' in dopts:
            kscalc.chkfile = dopts['chkfl']
        else:
            kscalc.chkfile = './WVFNS'

    if dopts['init'] == 'chkfl':
        kscalc.chkfile = dopts['chkfl']
        kscalc.init_guess = 'chkfile'
    elif dopts['init'] == 'HF DM':
        from pyscf import scf
        if dopts['restricted']:
            hfcalc = scf.RHF(mol)
        else:
            hfcalc = scf.UHF(mol)

        hfcalc.max_cycle = dopts['max_cycle']
        hfcalc.conv_tol = dopts['tol']
        if 'levelshift' in dopts:
            hfcalc.level_shift = dopts['levelshift']
        hfcalc.kernel()
        kscalc.init_guess = hfcalc.make_rdm1()

    kscalc.max_cycle = dopts['max_cycle']
    kscalc.conv_tol = dopts['tol']

    if dopts['xc'] != 'HF':

        kscalc.grids.level = dopts['gridsize']

        if dopts['xc_lib'] == 'XCFun':
            kscalc._numint.libxc = dft.xcfun
        elif dopts['xc_lib'] != 'LibXC':
            raise SystemExit('Unknown XC library '+ dopts['xc_lib'])
        kscalc.xc = dopts['xc']

    if 'levelshift' in dopts:
        kscalc.level_shift = dopts['levelshift']

    e0 = kscalc.kernel()

    odict = {
        'Etot': kscalc.e_tot, 'Converged': kscalc.converged,
    }
    fname = './pyscf_run.yaml'
    with open(fname,'w+') as tfl:
        for akey in odict:
            tfl.write('{:}: {:}\n'.format(akey,odict[akey]))

    np.save('./orbital_eigenvalues.npy',kscalc.mo_energy)
    np.save('./orbital_occupancies.npy',kscalc.mo_occ)
    np.save('./orbital_coefficients.npy',kscalc.mo_coeff)

    return e0

if __name__ == "__main__":

    #print(parse_inp('./sample.txt'))
    molscf()
