import numpy as np
from pyscf import gto,dft
from nrlmol_basis import get_dfo_basis

opt_type = {
    'float': ['tol','levelshift'],
    'int': ['gridsize','max_cycle','charge','2S','verbose'],
    'bool': ['symm','restricted'],
    'str': ['xyz_fl','basis','xc','xc_lib','ecp','logfile','init','chkfl']
}

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

    odict = {}

    tbas = get_dfo_basis('Ne',basis='dfo.gbs')

    mol = gto.M(atom='Ne 0.0 0.0 0.0', basis=tbas, symmetry=False, \
        charge=0, spin=0, output = './Ne.log', \
        verbose = 4)

    for grid in range(1,10,1):

        kscalc = dft.RKS(mol)

        kscalc.max_cycle = 500
        kscalc.conv_tol = 1.e-8
        kscalc.grids.level = grid

        kscalc.xc = 'LDA_X, LDA_C_PW'

        e0 = kscalc.kernel()

        odict[grid] = {
            'Etot': kscalc.e_tot, 'Converged': kscalc.converged,
        }

    fname = './pyscf_run.yaml'
    with open(fname,'w+') as tfl:
        for akey in odict:
            tfl.write('{:}:\n'.format(akey))
            for bkey in odict[akey]:
                tfl.write('  {:}: {:}\n'.format(bkey,odict[akey][bkey]))

    return e0

if __name__ == "__main__":

    #print(parse_inp('./sample.txt'))
    molscf()
