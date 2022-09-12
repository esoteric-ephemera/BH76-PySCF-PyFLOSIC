from pyscf.gto import load as gtoload

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

def get_nrlmol_bas(xyz_file,var='DFO'):

    if var == 'DFO':
        wfl = './dfo-nrlmol.dat'
    elif var == 'DFO+':
        wfl = './dfo+-nrlmol.dat'
    else:
        raise ValueError('Unknown basis set {:}; stopping.'.format(var))

    atd = get_ats_from_xyz(xyz_file)
    basd = {}
    for at in atd:
        basd[at] = gtoload(wfl,symb=at)
    return basd

if __name__ == "__main__":

    bd = get_nrlmol_bas('./BH76_RKT22.xyz',var='DFO')
    print(bd)
