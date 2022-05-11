import yaml
from os import sys,getcwd,listdir,path,system

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

"""
rdir is the directory where the files
         BH76_ref_energies.yaml and BH76_chg_2s.yaml
are located
"""
rdir = '/Users/aaronkaplan/Dropbox/phd.nosync/computation/BH76/'

def BH76_analysis(dir='./',edict={},fprefix=''):

    refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
        Loader=yaml.Loader)

    cmd = yaml.load(open(rdir + 'BH76_chg_2s.yaml','r'),\
        Loader=yaml.Loader)

    if len(edict.keys()) == 0:
        edict = get_en_dict_yaml(cmd,wdir=dir)

    rxd = {}
    rxd['SP'] = edict # single-point
    rxd['RX'] = {} # reactions

    ofl = dir + fprefix + 'BH76_individual.csv'
    ostr = 'System, Energy (Hartree)\n'
    for asys in edict:
        ostr += '{:}, {:}\n'.format(asys,edict[asys])

    with open(ofl,'w+') as lfl:
        lfl.write(ostr)

    md = 0.0
    mad = 0.0
    rmsd = 0.0
    amax = [-1e20,None]
    amin = [1e20,None]

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
        rmsd += aerr**2
        if aerr > amax[0]:
            amax[0] = aerr
            amax[1] = indx
        if aerr < amin[0]:
            amin[0] = aerr
            amin[1] = indx
        #amax = max(amax,aerr)
        #amin = min(amin,aerr)

        ostr += ('{:}, '*3 + '{:}\n').format(indx,tmp,refs[indx]['Ref'],err)

        nrx += 1

    md /= (1.0*nrx)
    mad /= (1.0*nrx)
    rmsd = max(0.0,rmsd/(1.0*nrx) - mad**2)**(0.5)
    ostr += '--,--,--,--\n'
    ostr += 'MD (kcal/mol), {:}, , \n'.format(md)
    ostr += 'MAD (kcal/mol), {:}, , \n'.format(mad)
    ostr += 'RMSD (kcal/mol), {:}, , \n'.format(rmsd)
    ostr += 'AMAX (kcal/mol), {:}, index, {:} \n'.format(*amax)
    ostr += 'AMIN (kcal/mol), {:}, index, {:} '.format(*amin)

    ofl = dir + fprefix + 'BH76_errors.csv'
    with open(ofl,'w+') as lfl:
        lfl.write(ostr)

    rxd['Stats'] = {'MD': md, 'MAD': mad, 'RMSD': rmsd,
        'AMAX': {'Value': amax[0], 'Index': amax[1]},
        'AMIN': {'Value': amin[0], 'Index': amin[1]}
    }

    ofl = dir + fprefix + 'BH76_total.yaml'
    with open(ofl,'w+') as lfl:
        yaml.dump(rxd,lfl,Dumper=yaml.Dumper)

    return


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
