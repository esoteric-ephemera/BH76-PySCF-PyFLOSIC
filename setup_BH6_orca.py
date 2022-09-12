import yaml
from os import path, system, sys
from subprocess import run as subrun

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

bh6_rxs = [45,46,61,62,63,64]
rdir = path.dirname(path.realpath(__file__)) + '/'

accept_opts = ['method', 'basis', '*xyz', '*xyzfile', '%scf', 'charge', '2S+1',
    'tol', 'grid', 'rest', '%pal']

def get_BH6_sys_list():
    rxd = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
        Loader=yaml.Loader)

    sysl = []
    for irx in bh6_rxs:
        td = rxd[irx]['Stoich']
        for asys in td:
            if asys not in sysl:
                sysl.append(asys)
    return sysl

def job_script(file_name='./runjob.sh'):

    tstr = '#!/bin/bash\n'
    tstr += '#PBS -l walltime=1:00:00:00\n'
    tstr += '#PBS -q normal\n'
    tstr += '#PBS -l nodes=1:ppn=20\n'
    tstr += '#PBS -N BH6\n'
    tstr += '#PBS -j oe\n'
    tstr += '#PBS -o oe.txt\n'
    tstr += '#PBS -m a\n\n'

    tstr += 'cd "$PBS_O_WORKDIR"\n'
    tstr += 'export PATH=/home/tuf53878/orca/openmpi-4.1.1/bin:$PATH\n'
    tstr += 'export LD_LIBRARY_PATH=/home/tuf53878/orca/openmpi-4.1.1/lib:$LD_LIBRARY_PATH\n'
    tstr += 'orca_exec=/home/tuf53878/orca/orca\n\n'

    tstr += 'for u in ./*/ ; do\n'
    tstr += '  cd $u\n'
    tstr += '  $orca_exec inp.txt >> out.txt\n'
    tstr += '  cd ..\n'
    tstr += 'done'


    with open(file_name,'w+') as tfl:
        tfl.write(tstr)

    return

def orca_inp_writer(opts={},file_name='./inp.txt'):
    def_opts = {
        'method': 'PBE', 'basis': 'def2-qzvp', '*xyzfile': 'struc.xyz',
        '%scf': {'MaxIter': 500}, 'charge': 0, '2S+1': 1, 'tol': 'StrongSCF',
        'grid': 'defgrid2', 'rest': 'RKS', '%pal': {'nprocs': 20}
    }

    for anopt in opts:
        def_opts[anopt] = opts[anopt]

    with open(file_name,'w+') as tfl:
        for anopt in def_opts:
            if anopt[0] == '%':
                tfl.write('{:}\n'.format(anopt))
                for bopt in def_opts[anopt]:
                    tfl.write('{:}  {:}\n'.format(bopt,def_opts[anopt][bopt]))
                tfl.write('end\n')
            elif anopt[0] == '*':
                if anopt == '*xyz':
                    tfl.write('{:}  {:}  {:}\n'.format(anopt,def_opts['charge'],def_opts['2S+1']))
                    for anat in def_opts[anopt]:
                        tfl.write('{:}  {:}  {:}  {:}\n'.format(anat,*def_opts[anopt][anat]))
                    tfl.write('*\n')
                elif anopt == '*xyzfile':
                    tfl.write('{:}  {:}  {:}  {:}\n'.format(anopt, def_opts['charge'], \
                        def_opts['2S+1'], def_opts['*xyzfile']))

            elif anopt in ['charge','2S+1']:
                continue
            else:
                tfl.write('! {:}\n'.format(def_opts[anopt]))

    return

def setup_BH6(uopts,makeclean=True):

    cmd = yaml.load(open(rdir + 'BH76_chg_2s.yaml','r'),\
        Loader=yaml.Loader)

    sys_l = get_BH6_sys_list()

    thstr = 'KS'
    if uopts['method'] == 'HF':
        thstr = 'HF'

    wdir = 'BH6_{:}/'.format(uopts['method'])
    geodir = rdir + 'BH76_geometries/'

    if makeclean:
        system('rm -rf ' + wdir)
    system('mkdir ' + wdir)

    job_script(file_name=wdir+'runjob.sh')

    for asys in sys_l:
        system('cp -r {:}/{:} {:}/'.format(geodir,asys,wdir))

        topts = { }
        for anopt in uopts:
            topts[anopt] = uopts[anopt]
        topts['charge'] = cmd[asys]['Charge']
        topts['2S+1'] = cmd[asys]['2S']+1
        if cmd[asys]['R']:
            topts['rest'] = 'R'+thstr
        else:
            topts['rest'] = 'U'+thstr

        orca_inp_writer(opts=topts,file_name=wdir + asys +'/inp.txt')

    return

def analysis(dir='./'):

    rxd = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
        Loader=yaml.Loader)

    sysl = get_BH6_sys_list()
    en_d = {
        'units': 'kcal/mol',
        'SP': {},
        'RX': {},
        'Stats': {'MD': 0.0, 'MAD': 0.0, 'VAR': 0.0, 'STD DEV': 0.0, 'RMSD': 0.0}
    }
    for asys in sysl:
        tmp_byte = subrun(['grep','FINAL SINGLE POINT ENERGY',asys+'/out.txt'], \
            capture_output=True)
        vals = tmp_byte.stdout.decode('utf-8').split()
        etot = float(vals[-1])
        tmp_byte = subrun(['grep','SCF CONVERGED',asys+'/out.txt'], \
            capture_output=True)
        tmpstr = tmp_byte.stdout.decode('utf-8').split()
        if len(tmpstr) == 0:
            print('WARNING, {:} not converged'.format(asys))
        en_d['SP'][asys] = etot*eH_to_kcalmol

    for jrx,irx in enumerate(bh6_rxs):
        en_d['RX'][jrx] = 0.0
        td = rxd[irx]['Stoich']
        for asys in td:
            en_d['RX'][jrx] += td[asys]*en_d['SP'][asys]
        terr = en_d['RX'][jrx] - rxd[irx]['Ref']

        en_d['Stats']['MD'] += terr
        en_d['Stats']['MAD'] += abs(terr)
        en_d['Stats']['VAR'] += terr**2

    fnsys = 1.*len(bh6_rxs)
    en_d['Stats']['MD'] /= fnsys
    en_d['Stats']['MAD'] /= fnsys
    en_d['Stats']['VAR'] /= fnsys
    en_d['Stats']['STD DEV'] = max(0.0, en_d['Stats']['VAR'] - en_d['Stats']['MD']**2)**(0.5)
    en_d['Stats']['RMSD'] = en_d['Stats']['VAR']**(0.5)

    with open('./BH6_total.yaml','w+') as tfl:
        yaml.dump(en_d,tfl,Dumper=yaml.Dumper)

    return

def parse_into_setup():

    method = None
    makeclean = True
    dfa_dir = None

    ioptd = {}

    if len(sys.argv) < 2:
        raise SystemExit('Need to specify at least the density functional!')
    elif sys.argv[1] == 'analysis':
        analysis()
        exit()
    else:
        for astr in sys.argv[1:]:
            key,val = astr.split('=')
            lkey = key.lower()
            if lkey == 'clean':
                if val.lower()[0] == 't':
                    makeclean = True
                elif val.lower()[0] == 'f':
                    makeclean = False
                else:
                    print('Keyword "clean" must be t/f!')
            elif key in accept_opts:
                ioptd[key] = val
            else:
                print('Skipping unknown keyword {:}'.format(key))

    setup_BH6(ioptd,makeclean=makeclean)

    return


if __name__ == "__main__":

    parse_into_setup()
