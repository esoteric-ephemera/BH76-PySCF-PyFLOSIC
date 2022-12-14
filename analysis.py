import yaml
from os import sys,getcwd,listdir,path,system

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

bh6_rxs = [45,46,61,62,63,64]

"""
rdir is the directory where the files
         BH76_ref_energies.yaml and BH76_chg_2s.yaml
are located
"""

rdir = path.dirname(path.realpath(__file__)) + '/'

def BH76_analysis(cdir='./',edict={},fprefix='',wrc=True):

    if cdir[-1] != '/':
        cdir += '/'

    refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
        Loader=yaml.Loader)

    cmd = yaml.load(open(rdir + 'BH76_chg_2s.yaml','r'),\
        Loader=yaml.Loader)

    if len(edict.keys()) == 0:
        edict = get_en_dict_yaml(cmd,wdir=cdir)

    rxd = {}
    rxd['SP'] = edict # single-point
    rxd['RX'] = {} # reactions

    ofl = cdir + fprefix + 'BH76_individual.csv'
    ostr = 'System, Energy (Hartree)\n'
    for asys in edict:
        ostr += '{:}, {:}\n'.format(asys,edict[asys])

    with open(ofl,'w+') as lfl:
        lfl.write(ostr)

    md = 0.0
    mad = 0.0
    var = 0.0
    mape = 0.0
    amax = [-1e20,None]
    amin = [1e20,None]

    bh6_md = 0.0
    bh6_mad = 0.0
    bh6_var = 0.0
    bh6_mape = 0.0

    ostr = 'Reaction index, Energy (kcal/mol), Ref (kcal/mol), Error (kcal/mol)\n'
    nrx = 0
    for indx in refs:

        rxd['RX'][indx] = {}

        tmp = 0.0
        td = refs[indx]['Stoich']
        for aspec in td:
            if aspec in edict:
                tmp += td[aspec]*edict[aspec]
            else:
                print('WARNING, no computed energy for system '+aspec)

        tmp *= eH_to_kcalmol
        err = tmp - refs[indx]['Ref']
        rxd['RX'][indx]['Energy'] = tmp
        rxd['RX'][indx]['Error'] = err
        aerr = abs(err)
        md += err
        mad += aerr
        var += aerr**2
        mape += 100*aerr/abs(refs[indx]['Ref'])

        if aerr > amax[0]:
            amax[0] = aerr
            amax[1] = indx
        if aerr < amin[0]:
            amin[0] = aerr
            amin[1] = indx
        if indx in bh6_rxs:
            bh6_md += err
            bh6_mad += aerr
            bh6_var += aerr**2
            bh6_mape += 100*aerr/abs(refs[indx]['Ref'])

        ostr += ('{:}, '*3 + '{:}\n').format(indx,tmp,refs[indx]['Ref'],err)

        nrx += 1

    fnrx = 1.0*nrx
    md /= fnrx
    mad /= fnrx
    var /= fnrx
    mape /= fnrx
    rmsd = var**(0.5)
    stddev = max(0.0,var - mad**2)**(0.5)
    ostr += '--,--,--,--\n'
    ostr += 'MD (kcal/mol), {:}, , \n'.format(md)
    ostr += 'MAD (kcal/mol), {:}, , \n'.format(mad)
    ostr += 'RMSD (kcal/mol), {:}, , \n'.format(rmsd)
    ostr += 'VAR (kcal/mol), {:}, , \n'.format(var)
    ostr += 'STDDEV (kcal/mol), {:}, , \n'.format(stddev)
    ostr += 'MAPE (%), {:}, , \n'.format(mape)
    ostr += 'AMAX (kcal/mol), {:}, index, {:} \n'.format(*amax)
    ostr += 'AMIN (kcal/mol), {:}, index, {:} \n'.format(*amin)

    bh6_md /= 6.0
    bh6_mad /= 6.0
    bh6_var /= 6.0
    bh6_mape /= 6.0
    bh6_rmsd = bh6_var**(0.5)
    bh6_stddev = max(0.0,bh6_var - bh6_mad**2)**(0.5)
    ostr += '\nBH6 MD (kcal/mol), {:}, , \n'.format(bh6_md)
    ostr += 'BH6 MAD (kcal/mol), {:}, , \n'.format(bh6_mad)
    ostr += 'BH6 RMSD (kcal/mol), {:}, , \n'.format(bh6_rmsd)
    ostr += 'BH6 VAR (kcal/mol), {:}, , \n'.format(bh6_var)
    ostr += 'BH6 STDDEV (kcal/mol), {:}, , \n'.format(bh6_stddev)
    ostr += 'BH6 MAPE (%), {:}, , \n'.format(bh6_mape)

    ofl = cdir + fprefix + 'BH76_errors.csv'
    with open(ofl,'w+') as lfl:
        lfl.write(ostr)

    rxd['Stats'] = {'MD': md, 'MAD': mad, 'RMSD': rmsd,
        'VAR': var, 'STDDEV': stddev, 'MAPE': mape,
        'AMAX': {'Value': amax[0], 'Index': amax[1]},
        'AMIN': {'Value': amin[0], 'Index': amin[1]}
    }

    rxd['BH6 Stats'] = {'MD': bh6_md, 'MAD': bh6_mad, 'RMSD': bh6_rmsd,
        'VAR': bh6_var, 'STDDEV': bh6_stddev, 'MAPE': bh6_mape
    }

    if wrc:
        # also do BH76RC - no extra calculations needed

        refs = yaml.load(open(rdir + 'BH76RC_ref_energies.yaml','r'),\
            Loader=yaml.Loader)

        rxd['RC'] = {} # reactions

        md = 0.0
        mad = 0.0
        var = 0.0
        mape = 0.0
        amax = [-1e20,None]
        amin = [1e20,None]

        ostr = 'Reaction index, Energy (kcal/mol), Ref (kcal/mol), Error (kcal/mol)\n'
        nrx = 0
        for indx in refs:

            rxd['RC'][indx] = {}

            tmp = 0.0
            td = refs[indx]['Stoich']
            for aspec in td:
                if aspec in edict:
                    tmp += td[aspec]*edict[aspec]
                else:
                    print('WARNING, no computed energy for system '+aspec)
            tmp *= eH_to_kcalmol
            err = tmp - refs[indx]['Ref']
            rxd['RC'][indx]['Energy'] = tmp
            rxd['RC'][indx]['Error'] = err
            aerr = abs(err)
            md += err
            mad += aerr
            var += aerr**2
            mape += 100*aerr/abs(refs[indx]['Ref'])
            if aerr > amax[0]:
                amax[0] = aerr
                amax[1] = indx
            if aerr < amin[0]:
                amin[0] = aerr
                amin[1] = indx
            if indx in bh6_rxs:
                bh6_md += err
                bh6_mad += aerr
                bh6_var += aerr**2

            ostr += ('{:}, '*3 + '{:}\n').format(indx,tmp,refs[indx]['Ref'],err)

            nrx += 1

        fnrx = 1.0*nrx
        md /= fnrx
        mad /= fnrx
        var /= fnrx
        mape /= fnrx
        rmsd = var**(0.5)
        stddev = max(0.0,var - mad**2)**(0.5)
        ostr += '--,--,--,--\n'
        ostr += 'MD (kcal/mol), {:}, , \n'.format(md)
        ostr += 'MAD (kcal/mol), {:}, , \n'.format(mad)
        ostr += 'RMSD (kcal/mol), {:}, , \n'.format(rmsd)
        ostr += 'VAR (kcal/mol), {:}, , \n'.format(var)
        ostr += 'STDDEV (kcal/mol), {:}, , \n'.format(stddev)
        ostr += 'MAPE (%), {:}, , \n'.format(mape)
        ostr += 'AMAX (kcal/mol), {:}, index, {:} \n'.format(*amax)
        ostr += 'AMIN (kcal/mol), {:}, index, {:} \n'.format(*amin)

        ofl = cdir + fprefix + 'BH76RC_errors.csv'
        with open(ofl,'w+') as lfl:
            lfl.write(ostr)

        rxd['RC Stats'] = {'MD': md, 'MAD': mad, 'RMSD': rmsd,
            'VAR': var, 'STDDEV': stddev, 'MAPE': mape,
            'AMAX': {'Value': amax[0], 'Index': amax[1]},
            'AMIN': {'Value': amin[0], 'Index': amin[1]}
        }

    ofl = cdir + fprefix + 'BH76_total.yaml'
    with open(ofl,'w+') as lfl:
        yaml.dump(rxd,lfl,Dumper=yaml.Dumper)

    return rxd


def get_en_dict_yaml(refed,wdir='./'):

    if path.isfile(wdir + '/non_conv_list.txt'):
        system('rm ' + wdir + '/non_conv_list.txt')

    non_conv = 0
    nc_list = []
    en_dict = {}

    for adir in refed:#listdir(wdir):

        if not path.isdir(wdir+'/'+adir) or adir[0]=='.':
            continue


        tmpf = wdir+'/'+adir + '/pyscf_run.yaml'
        if path.isfile(tmpf):
            tmpd = yaml.load(open(tmpf,'r'),Loader=yaml.Loader)
        else:
            tmpd = {'Etot': 0.0, 'Converged': False}

        en_dict[adir] = tmpd['Etot']
        if not tmpd['Converged']:
            non_conv += 1
            nc_list.append(adir)

    if non_conv > 0:
        with open(wdir + '/non_conv_list.txt','w+') as tfl:
            tfl.write('Total = {:}\n'.format(non_conv))
            for tmp in nc_list:
                tfl.write('{:}\n'.format(tmp))

    return en_dict

if __name__ == "__main__":

    BH76_analysis()
