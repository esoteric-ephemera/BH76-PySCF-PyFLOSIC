import yaml
from os import path,system
from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

#yaml.Dumper.ignore_aliases = lambda *args : True

conversion() # ensures FLOSIC errors are recomputed
no_flo_dat = ['LCwPBE','BLYP','B3LYP', 'M06-L','MN15-L']

rdir = path.dirname(path.realpath(__file__)) + '/'

refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'), Loader=yaml.Loader)
cmd = yaml.load(open(rdir + 'BH76_chg_2s.yaml','r'), Loader=yaml.Loader)

def get_wgts():
    rpf = {}
    rpb = {}
    ts = {}
    for ibh, abh in enumerate(refs):
        for asys in refs[abh]['Stoich']:
            if refs[abh]['Stoich'][asys] > 0:
                wdct = ts
            else:
                if ibh%2 == 0:
                    wdct = rpf
                else:
                    wdct = rpb

            if asys in wdct:
                wdct[asys] += 1
            else:
                wdct[asys] = 1

    for wdct in [rpf, rpb]:
        for akey in wdct:
            wdct[akey] /= 38.

    for akey in ts:
        ts[akey] /= 76.

    return ts, rpf, rpb

ts_d, rpf_d, rpb_d = get_wgts()

out_dir = './dc_dft_analysis/'
if not path.isdir(out_dir):
    system('mkdir '+out_dir)

scf_dir = './results_aug-cc-pvqz/'
nrlmol_basis_dir = './results_NRLMOL_cart/'
refs = {
    #'SCAN@HF': './results_aug-cc-pvqz/',
    'SCAN-FLOSIC': './FLOSIC/',
    'LCwPBE': './results_aug-cc-pvqz/',
    'S50X': './scan_hybs/S50X/'
    #'S25X': './scan_hybs/S25X/',
}
Nref = len(refs.keys())

refd = {}

dfas = ['LSDA','PBE','BLYP','SCAN','r2SCAN','M06-L','MN15-L','B3LYP','LCwPBE','SCAN@HF']
metrics = ['MFE','MAFE','MDE','MADE']

for aref in refs:
    if aref in ['S25X','S50X']:
        frac = aref[1:-1]
        tdir = refs[aref]+'SCAN_{:}_EXX_BH76/'.format(frac)
    elif aref in ['SCAN@HF','SCAN-FLOSIC']:
        tdir = refs[aref] + 'SCAN/' + aref +'_BH76/'
    else:
        tdir = refs[aref] + aref + '/' + aref +'_BH76/'

    tfl = tdir + 'BH76_total.yaml'
    if aref != 'SCAN-FLOSIC': # conversion() takes care of SCAN-FLOSIC
        BH76_analysis(cdir=tdir)
    refd[aref] = yaml.load(open(tfl,'r'), Loader=yaml.Loader)#['SP']

tstr = {'TS': '', 'RPF': '', 'RPB': '', 'BH': '', 'RC': ''}

for akey in tstr:
    tstr[akey] += akey

for amet in metrics:
    for aref in refs:
        for akey in tstr:
            tstr[akey] += ', {:}'.format(aref)
for akey in tstr:
    tstr[akey] += '\n'

for dfa in dfas:

    for akey in tstr:
        tstr[akey] += '{:}'.format(dfa)

    wd = {}
    for amet in metrics:
        wd[amet] = {}
        for aref in refs:
            wd[amet][aref] = {
                'TS': 0.0, 'RPF': 0.0, 'RPB': 0.0, 'BH': 0.0, 'RC': 0.0
            }

    for aref in refs:

        if dfa == aref:
            continue
        if (aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat):
            continue

        if aref == 'SCAN-FLOSIC':

            if dfa == 'SCAN@HF':
                tdir_scf = nrlmol_basis_dir+'SCAN/SCAN@HF_BH76/'
                tdir_nscf = './FLOSIC/SCAN@{:}_BH76/'.format(aref)
            else:
                tdir_scf = './FLOSIC/{:}/{:}_BH76/'.format(dfa,dfa)
                tdir_nscf = './FLOSIC/{:}/{:}@{:}_BH76/'.format(dfa,dfa,aref)

            td_scf = yaml.load(open(tdir_scf + 'BH76_total.yaml','r'), Loader=yaml.Loader)
            td_nscf = yaml.load(open(tdir_nscf+'BH76_total.yaml','r'), Loader=yaml.Loader)

        else:

            if dfa == 'SCAN@HF':
                tdir = scf_dir+'SCAN/SCAN@HF_BH76/'
            else:
                tdir = scf_dir+'{:}/{:}_BH76/'.format(dfa,dfa)

            tfl = tdir + 'BH76_total.yaml'
            BH76_analysis(cdir=tdir)
            td_scf = yaml.load(open(tfl,'r'), Loader=yaml.Loader)

            tref = aref
            if aref == 'SCAN@HF':
                tref = 'HF'

            if aref in ['S25X','S50X']:
                if dfa == 'SCAN@HF':
                    tdir = refs[aref]+'SCAN@{:}_BH76/'.format(tref)
                else:
                    tdir = refs[aref]+'{:}@{:}_BH76/'.format(dfa,tref)
            else:
                if dfa == 'SCAN@HF':
                    tdir = refs[aref]+'SCAN/SCAN@{:}_BH76/'.format(tref)
                else:
                    tdir = refs[aref]+'{:}/{:}@{:}_BH76/'.format(dfa,dfa,tref)
            #print(aref,dfa,tdir)
            tfl = tdir + 'BH76_total.yaml'
            BH76_analysis(cdir=tdir)
            td_nscf = yaml.load(open(tfl,'r'), Loader=yaml.Loader)

        for idict,adict in enumerate([ts_d,rpf_d,rpb_d]):

            if idict == 0:
                didx = 'SP'
                idx = 'TS'
            elif idict == 1:
                didx = 'SP'
                idx = 'RPF'
            elif idict == 2:
                didx = 'SP'
                idx = 'RPB'

            for asys in adict:

                tfe = (td_nscf[didx][asys] - refd[aref][didx][asys])*eH_to_kcalmol
                tde = (td_scf[didx][asys] - td_nscf[didx][asys])*eH_to_kcalmol

                wd['MFE'][aref][idx] += tfe*adict[asys]
                wd['MAFE'][aref][idx] += abs(tfe)*adict[asys]
                wd['MDE'][aref][idx] += tde*adict[asys]
                wd['MADE'][aref][idx] += abs(tde)*adict[asys]

        for idict,adict in enumerate([td_scf['RX'],td_scf['RC']]):

            if idict == 0:
                didx = 'RX'
                idx = 'BH'
            elif idict == 1:
                didx = 'RC'
                idx = 'RC'

            for asys in td_scf[didx]:
                tfe = (td_nscf[didx][asys]['Energy'] - refd[aref][didx][asys]['Energy'])
                tde = (td_scf[didx][asys]['Energy'] - td_nscf[didx][asys]['Energy'])

                wd['MFE'][aref][idx] += tfe
                wd['MAFE'][aref][idx] += abs(tfe)
                wd['MDE'][aref][idx] += tde
                wd['MADE'][aref][idx] += abs(tde)

        for amet in metrics:
            wd[amet][aref]['BH'] /= (1.0*len(td_scf['RX'].keys()))
            wd[amet][aref]['RC'] /= (1.0*len(td_scf['RC'].keys()))

    ofl = out_dir + dfa + '_BH76_dc_dft.yaml'
    with open(ofl,'w+') as tfl:
        yaml.dump(wd,tfl,Dumper=yaml.Dumper)

    for amet in metrics:
        for aref in refs:
            for akey in tstr:
                tstr[akey] += ', {:}'.format(wd[amet][aref][akey])
    for akey in tstr:
        tstr[akey] += '\n'

with open(out_dir+'dc_dft_analysis.csv','w+') as tfl:
    for akey in tstr:
        tfl.write(tstr[akey])

lbld = {'BH': 'BH76', 'RC': 'BH76RC', 'RPF': 'Reactants and products',
    'RPB': 'Reactants and products',
    'TS': 'Transition states'
}

for i in range(2):
    if i == 0:
        nstr = {'TS': '', 'RPF': '', 'RPB': ''}
        tchar = 'SP'
    elif i == 1:
        nstr = {'BH': '', 'RC': ''}
        tchar = 'delta_E'

    for akey in nstr:
        if akey == 'RPF':
            nstr[akey] += '\\textit{Forwards}'
        elif akey == 'RPB':
            nstr[akey] += '\\textit{Reverse}'

        for aref in refd:
            if aref == 'LCwPBE':
                aref_str = '\\lcwpbe{}'
            else:
                aref_str = aref
            nstr[akey] += ' & \\multicolumn{2}{c}{'+'{:}'.format(aref_str) + '}'
        nstr[akey] += ' \\\\ \n'

        nstr[akey] += '\\textit{'+lbld[akey]+'}'

        for aref in refd:
            nstr[akey] += ' & MFE & MDE'
        nstr[akey] += ' \\\\ \\hline \n'

    for idfa,dfa in enumerate(dfas):
        ofl = out_dir + dfa + '_BH76_dc_dft.yaml'
        wd = yaml.load(open(ofl,'r'), Loader=yaml.Loader)

        for akey in nstr:
            if dfa == 'r2SCAN':
                dfa_str = '\\rrscan{}'
            elif dfa == 'LCwPBE':
                dfa_str = '\\lcwpbe{}'
            else:
                dfa_str = dfa
            nstr[akey] += dfa_str

        for aref in refd:
            for akey in nstr:
                #pe_rat = 100*wd['MADE'][aref][akey]/max(1.e-15,wd['MAFE'][aref][akey])
                if (dfa == aref) or \
                    ( (aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat)):
                    nstr[akey] += ' &  & '
                else:
                    nstr[akey] += ' & {:.2f} & {:.2f}'.format(wd['MFE'][aref][akey],wd['MDE'][aref][akey])

        lchar = ''
        if idfa == len(dfas)-1:
            lchar = '\\hline'

        for akey in nstr:
            nstr[akey] += ' \\\\ '+lchar+' \n'

    with open(out_dir+'mde_mfe_{:}.tex'.format(tchar),'w+') as tfl:
        for akey in nstr:
            tfl.write(nstr[akey])
