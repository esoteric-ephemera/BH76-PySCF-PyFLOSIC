from setup import setup
from glob import glob
from os import system, path,listdir

def get_wvfns(dfa,basis='aug-cc-pvqz',bdir='./'):
    tdir = '{:}{:}_BH76/{:}_WVFN_{:}/'.format(bdir,dfa,dfa,basis)
    if not path.isdir(tdir):
        system('mkdir '+tdir)

    for asys in listdir('{:}{:}_BH76/'.format(bdir,dfa)):
        tsys = '{:}{:}_BH76/{:}'.format(bdir,dfa,asys)
        if path.isdir(tsys) and '_WVFN' not in asys:
            system('cp {:}/WVFNS {:}/{:}_wvfn'.format(tsys,tdir,asys))
    return

def setup_hyb_calcs(hdir,sl_name,basis='aug-cc-pvqz'):
    dfas = ['LSDA', 'PBE', 'BLYP' , 'SCAN', 'r2SCAN', 'B3LYP' ,'LCwPBE']

    tdir_l = glob(hdir+'*_EXX_BH76/')

    dirs = [sl_name]
    for tdir in tdir_l:
        tmp = (tdir.split('/')[-2])[:-5]
        dirs.append(tmp)

    for tdir in dirs:

        if tdir == sl_name:
            frac = '00'
        else:
            frac = tdir.split('_EXX')[0]
            frac = frac.split('{:}_'.format(sl_name))[1]

        get_wvfns(tdir,basis=basis,bdir=hdir)

        wdir = '{:}/{:}{:}X/'.format(hdir,sl_name,frac)
        system('mkdir ' + wdir)
        system('mv {:}{:}_BH76 {:}/{:}_BH76'.format(hdir,tdir,wdir,tdir))

        for dfa in dfas:

            if dfa == sl_name and frac == '00':
                continue

            dfa_dir = wdir
            if dfa == 'HYB_GGA_XC_LC_WPBE':
                dfa_dir += 'LCwPBE'
            else:
                dfa_dir += dfa
            dfa_dir += '@{:}{:}X'.format(sl_name,frac)

            at_dfa = wdir+ '{:}_BH76/{:}_WVFN_{:}'.format(tdir,tdir,basis)

            xc_str = dfa
            if frac == '00':
                xc_str += '@{:}'.format(sl_name)
            else:
                xc_str += '@{:}{:}X'.format(sl_name,frac)

            setup(xc_str,startclean=True,xclib='LibXC', dfa_dir=dfa_dir, \
                wfdir=at_dfa, inp_opts={'basis': basis})

    return

def add_new_calcs(dfas,hdir,sl_name,skip_dirs=[],basis='aug-cc-pvqz'):
    # use to add new @hybrid calculations to existing directories

    dirs = glob(hdir+'*X')
    adjdirs = [adir.split('/')[-1] for adir in dirs]

    for idir, wdir in enumerate(dirs):

        if adjdirs[idir] in skip_dirs:
            continue

        frac = (wdir[:-1]).split(sl_name)[-1]

        for dfa in dfas:

            if dfa == sl_name and frac == '00':
                continue

            dfa_dir = wdir + '/'
            if dfa == 'HYB_GGA_XC_LC_WPBE':
                dfa_dir += 'LCwPBE'
            else:
                dfa_dir += dfa
            dfa_dir += '@{:}{:}X'.format(sl_name,frac)

            if frac == '00':
                at_dfa = wdir + '/{:}_BH76/{:}_WVFN_{:}'.format(sl_name,sl_name,basis)
            else:
                at_dfa = wdir + '/{:}_{:}_EXX_BH76/{:}_{:}_EXX_WVFN_{:}'.format(\
                    sl_name,frac,sl_name,frac,basis)

            xc_str = dfa
            if frac == '00':
                xc_str += '@{:}'.format(sl_name)
            else:
                xc_str += '@{:}{:}X'.format(sl_name,frac)

            setup(xc_str,startclean=True,xclib='LibXC', dfa_dir=dfa_dir, \
                wfdir=at_dfa, inp_opts={'basis': basis})

    return


if __name__ == "__main__":

    #get_wvfns('LCwPBE',basis='aug-cc-pvqz')

    hdir = './r2SCAN_hybrids/'
    sl_name = 'r2SCAN'
    #setup_hyb_calcs(hdir,sl_name,basis='aug-cc-pvqz')
    add_new_calcs(['BLYP','B3LYP'],hdir,sl_name,skip_dirs=['r2SCAN00X'],basis='aug-cc-pvqz')
