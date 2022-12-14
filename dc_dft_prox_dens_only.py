import yaml
from os import path,system
from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

#yaml.Dumper.ignore_aliases = lambda *args : True

def rep_rules(instr):
    rules = {'r2SCAN': '\\rrscan{}','LCwPBE': '\\lcwpbe{}','S50X': 'SX-0.5'}
    outstr = instr
    for arule in rules:
        outstr = outstr.replace(arule,rules[arule])
    return outstr


conversion() # ensures FLOSIC errors are recomputed
no_flo_dat = ['LCwPBE','BLYP','B3LYP', 'M06-L','MN15-L']

rdir = path.dirname(path.realpath(__file__)) + '/'

refd = {}
for aset in ['BH76','BH76RC']:
    refd[aset] = {}
    td = yaml.load(open(rdir + aset + '_ref_energies.yaml','r'), \
        Loader=yaml.Loader)
    for arx in td:
        refd[aset][arx] = td[arx]['Ref']

out_dir = './dc_dft_prox_dens_only/'
if not path.isdir(out_dir):
    system('mkdir '+out_dir)

scf_dir = './results_aug-cc-pvqz/'
nrlmol_basis_dir = './results_NRLMOL_cart/'
refs = {
    #'SCAN@HF': './results_aug-cc-pvqz/',
    'LCwPBE': './results_aug-cc-pvqz/',
    'S50X': './scan_hybs/S50X/',
    'SCAN-FLOSIC': './FLOSIC/'
    #'S25X': './scan_hybs/S25X/',
}
Nref = len(refs.keys())

dfas = ['LSDA','PBE','BLYP','SCAN','r2SCAN','M06-L','MN15-L','B3LYP','LCwPBE','LSDA@HF','PBE@HF','SCAN@HF','r2SCAN@HF']
metrics = ['MFE','MAFE','MDE','MADE']

tstr = {'BH': '', 'RC': ''}
lbld = {'BH': 'BH76', 'RC': 'BH76RC'}

for akey in tstr:
    tstr[akey] += '{:} index'.format(akey)

for amet in metrics:
    for aref in refs:
        for akey in tstr:
            tstr[akey] += ', {:}'.format(aref)
for akey in tstr:
    tstr[akey] += '\n'

for dfa in dfas:

    kbstr = {}
    for aset in refd:
        kbstr[aset] = ['' for isys in range(len(refd[aset].keys())+1)]
        kbstr[aset][0] += '{:} index'.format(akey)
        for isys in range(len(refd[aset].keys())):
            kbstr[aset][isys+1] = '{:}'.format(isys+1)

    for akey in tstr:
        tstr[akey] += '{:}'.format(dfa)

    wd = {}
    for amet in metrics:
        wd[amet] = {}
        for aref in refs:
            wd[amet][aref] = {'BH': 0.0, 'RC': 0.0}

    dsplit = dfa.split('@')
    tdfa = dsplit[0]
    if len(dsplit) > 1:
        at_dfa = dsplit[1]

    tdir = scf_dir+'{:}/{:}_BH76/'.format(tdfa,dfa)

    tfl = tdir + 'BH76_total.yaml'
    td_scf_default = BH76_analysis(cdir=tdir)

    for aref in refs:

        if dfa == aref:
            continue
        if (aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat):
            continue
        for aset in kbstr:
            kbstr[aset][0] += '\t{:} FE\t{:} DE'.format(aref,aref)

        if aref == 'SCAN-FLOSIC':

            if '@' in dfa:
                tdir_scf = nrlmol_basis_dir+'{:}/{:}_BH76/'.format(tdfa,dfa)
                td_scf = BH76_analysis(cdir=tdir_scf)
            else:
                tdir_scf = './FLOSIC/{:}/{:}_BH76/'.format(dfa,dfa)
                td_scf = yaml.load(open(tdir_scf + 'BH76_total.yaml','r'),\
                    Loader=yaml.Loader)

            tdir_nscf = './FLOSIC/{:}/{:}@{:}_BH76/'.format(tdfa,tdfa,aref)
            td_nscf = yaml.load(open(tdir_nscf+'BH76_total.yaml','r'),\
                Loader=yaml.Loader)

        else:

            td_scf = td_scf_default.copy()

            tref = aref
            if aref == 'SCAN@HF':
                tref = 'HF'

            if aref in ['S25X','S50X']:
                tdir = refs[aref]+'{:}@{:}_BH76/'.format(tdfa,tref)
            else:
                tdir = refs[aref]+'{:}/{:}@{:}_BH76/'.format(tdfa,tdfa,tref)

            #print(aref,dfa,tdir)
            tfl = tdir + 'BH76_total.yaml'
            BH76_analysis(cdir=tdir)
            td_nscf = yaml.load(open(tfl,'r'), Loader=yaml.Loader)

        for idict in range(2):

            if idict == 0:
                didx = 'RX'
                idx = 'BH'
            elif idict == 1:
                didx = 'RC'
                idx = 'RC'

            for isys, asys in enumerate(td_scf[didx]):
                tfe = (td_nscf[didx][asys]['Energy'] - refd[lbld[idx]][asys])
                tde = (td_scf[didx][asys]['Energy'] - td_nscf[didx][asys]['Energy'])

                wd['MFE'][aref][idx] += tfe
                wd['MAFE'][aref][idx] += abs(tfe)
                wd['MDE'][aref][idx] += tde
                wd['MADE'][aref][idx] += abs(tde)

                if idict == 0:
                    aset = 'BH76'
                elif idict == 1:
                    aset = 'BH76RC'
                kbstr[aset][isys+1] += '\t{:.2f}\t{:.2f}'.format(tfe,tde)

        for amet in metrics:
            wd[amet][aref]['BH'] /= (1.0*len(td_scf['RX'].keys()))
            wd[amet][aref]['RC'] /= (1.0*len(td_scf['RC'].keys()))

    ofl = out_dir + dfa + '_BH76_dc_dft.yaml'
    with open(ofl,'w+') as tfl:
        yaml.dump(wd,tfl,Dumper=yaml.Dumper)

    if dfa in ['LSDA','PBE','SCAN']:
        for aset in kbstr:
            ofl = out_dir + dfa + '_{:}_dc_dft.tsv'.format(aset)
            with open(ofl,'w+') as tfl:
                for i in range(len(kbstr[aset])):
                    tfl.write(kbstr[aset][i] + '\n')

    for amet in metrics:
        for aref in refs:
            for akey in tstr:
                tstr[akey] += ', {:}'.format(wd[amet][aref][akey])
    for akey in tstr:
        tstr[akey] += '\n'

with open(out_dir+'dc_dft_analysis.csv','w+') as tfl:
    for akey in tstr:
        tfl.write(tstr[akey])

nstr = {'BH': '', 'RC': ''}
tchar = 'delta_E'

for akey in nstr:

    for aref in refs:
        if aref == 'LCwPBE':
            aref_str = '\\lcwpbe{}'
        else:
            aref_str = aref
        nstr[akey] += ' & \\multicolumn{2}{c}{'+'{:}'.format(aref_str) + '}'
    nstr[akey] += ' \\\\ \n'

    nstr[akey] += '\\textit{'+lbld[akey]+'}'

    for aref in refs:
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

    for aref in refs:
        for akey in nstr:
            #pe_rat = 100*wd['MADE'][aref][akey]/max(1.e-15,wd['MAFE'][aref][akey])
            if (dfa == aref) or \
                ( (aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat)):
                nstr[akey] += ' &  & '
            else:
                nstr[akey] += ' & {:.2f} & {:.2f}'\
                    .format(wd['MFE'][aref][akey],wd['MDE'][aref][akey])

    lchar = ''
    if idfa == len(dfas)-1:
        lchar = '\\hline'

    for akey in nstr:
        nstr[akey] += ' \\\\ '+lchar+' \n'

with open(out_dir+'mde_mfe_{:}.tex'.format(tchar),'w+') as tfl:
    for akey in nstr:
        tfl.write(rep_rules(nstr[akey]))
