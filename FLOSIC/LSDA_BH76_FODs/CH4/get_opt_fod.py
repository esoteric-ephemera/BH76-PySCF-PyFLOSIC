import numpy as np
import ase.io as aio
from nrlmol_basis import get_dfo_basis
from pycom import pycom_guess
from ase_pyflosic_optimizer import flosic_optimize
from flosic_os import xyz_to_nuclei_fod

opt_type = {
    'float': ['tol'],
    'int': ['gridsize','charge','2S','verbose'],
    'bool': ['symm'],
    'str': ['xyz_fl','basis','xc','ecp','logfile']
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
    return atds

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

def xyz_to_FRMORB(flnm,outfl='FRMORB'):

    [geo,nuclei,fod1,fod2,included] = xyz_to_nuclei_fod(aio.read(flnm))
    N_up = len(fod1)
    N_dn = len(fod2)
    fod_up = fod1.positions
    fod_dn = fod2.positions

    with open(outfl,'w+') as tfl:
        tfl.write('  {:} {:}\n'.format(N_up,N_dn))
        for i in range(N_up):
            tfl.write('{:.6f} {:.6f} {:.6f}\n'.format(*fod_up[i]))
        for i in range(N_dn):
            tfl.write('{:.6f} {:.6f} {:.6f}\n'.format(*fod_dn[i]))

    return

def fod_opt():

    dopts = {
        'gridsize': 5, 'basis': 'DFO', 'symm': False,
        'tol': 1.e-7, 'charge': 0, '2S': 0,
        'xc': 'LDA, PW', 'verbose': 4, 'ecp' : {}, 'opt': 'FIRE',
        'max_inc' : 0.1, 'label': 'OPT_FRMORB', 'fmax' : 1.e-4,
        'max_step' : 1000, 'max_cycle': 500
        }

    uopts = parse_inp('./inp.txt')

    for akey in uopts:
        dopts[akey] = uopts[akey]

    if dopts['basis'] in ['DFO','DFO+','UNC_DFO','UNC_DFO+']:
        tdict = get_ats_from_xyz(flnm)
        tstr = ''
        for at in tdict:
            tstr += at
        tbas = get_dfo_basis(tstr,basis=dopts['basis'].lower()+'.gbs')
    else:
        tbas = dopts['basis']

    pycom_guess(ase_nuclei=aio.read(dopts['xyz_fl']), charge = dopts['charge'], \
        spin = dopts['2S'], xc = dopts['xc'], basis = tbas, method = 'fb', \
        grid = dopts['gridsize'], calc ='UKS', symmetry = dopts['symm'],\
        verbose = dopts['verbose'])

    print('FODs found from PyCOM')
    xyz_to_FRMORB('FB_GUESS_COM.xyz',outfl='FRMORB_COM')

    fod_opt = flosic_optimize(mode='flosic-scf', atoms=aio.read('./FB_GUESS_COM.xyz'), \
        charge = dopts['charge'], spin = dopts['2S'], xc = dopts['xc'], \
        basis = tbas, ecp = dopts['ecp'], opt = dopts['opt'], \
        maxstep = dopts['max_inc'], label = dopts['label'], fmax = dopts['fmax'], \
        steps = dopts['max_step'], max_cycle = dopts['max_cycle'], conv_tol = dopts['tol'], \
        grid = dopts['gridsize'], ghost=False, use_newton=False, use_chk=False, \
        verbose=0, debug=False, efield=None, l_ij=None, ods=None, force_consistent=False, \
        fopt='force', fix_fods=False, ham_sic='HOOOV', vsic_every=1, n_rad=None, \
        n_ang=None, prune='no')

    print('FOD optimization complete')
    aio.write('FOD_OPT.xyz',fod_opt,format='xyz')
    xyz_to_FRMORB('FOD_OPT.xyz',outfl='FRMORB_OPT')

    return

if __name__ == "__main__":

    fod_opt()
