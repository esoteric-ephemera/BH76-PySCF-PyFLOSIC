import yaml
from os import path,system
import numpy as np

from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion
#from index_labels import indx_to_name

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

#yaml.Dumper.ignore_aliases = lambda *args : True

def str_rep_rules(instr):
    rules = {
        'r2SCAN': '\\rrscan{}',
        'LCwPBE': '\\lcwpbe{}',
        'S50X': 'SX-0.5'
    }
    outstr = instr
    for arule in rules:
        outstr = outstr.replace(arule,rules[arule])
    return outstr

# threshold for reliabilty of results
cutoff = 2.

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

fail_d = {}
for aset in refd:
    fail_d[aset] = [0 for i in range(len(refd[aset].keys()))]

out_dir = './dc_dft_cut/'
if not path.isdir(out_dir):
    system('mkdir '+out_dir)

scf_dir = './results_aug-cc-pvqz/'
nrlmol_basis_dir = './results_NRLMOL_cart/'
refs = {
    'LCwPBE': './results_aug-cc-pvqz/',
    'S50X': './scan_hybs/S50X/',
    'SCAN-FLOSIC': './FLOSIC/'
}
Nref = len(refs.keys())

dfas = ['LSDA','PBE','BLYP','SCAN','r2SCAN','M06-L','MN15-L','B3LYP','LCwPBE','LSDA@HF','PBE@HF','SCAN@HF','r2SCAN@HF']
#dfas = ['LSDA','LSDA@HF','PBE','PBE@HF','SCAN','SCAN@HF','r2SCAN','r2SCAN@HF']

metrics = ['MFE','MAFE','MDE','MADE']

tstr = {'TS': '', 'RPF': '', 'RPB': '', 'BH': '', 'RC': ''}

tex_str = {'BH76': '', 'BH76RC': ''}
for aset in tex_str:
    tex_str[aset] += ' & \\multicolumn{7}{c}{' + aset + '} \\\\ \n'
    tex_str[aset] += 'DFA & $\\overline{\\Delta E_\\mathrm{F}}$' \
        + ' & $\\overline{\\Delta E_\\mathrm{D}}$ ' \
        + ' & $\\delta \\widetilde{E}_\\mathrm{F}$ & $\\delta E_\\mathrm{D}$ ' \
        + ' & $N_\mathrm{cut}^{-1} \\sum_{j=1}^\\mathrm{N_\\mathrm{cut}} \\delta \\widetilde{E}_\\mathrm{F}^{(j)}$' \
        + ' & $N_\mathrm{cut}^{-1} \\sum_{j=1}^\\mathrm{N_\\mathrm{cut}} \\delta E_\\mathrm{D}^{(j)}$'\
        + '& No. passed \\\\ \hline \n'

for akey in tstr:
    tstr[akey] += akey

for amet in metrics:
    for aref in refs:
        for akey in tstr:
            tstr[akey] += ', {:}'.format(aref)
for akey in tstr:
    tstr[akey] += '\n'

#names = {'BH76': {}, 'BH76RC': {}}
#names['BH76'], names['BH76RC'] = indx_to_name()

for idfa, dfa in enumerate(dfas):

    failed = {'BH76': [], 'BH76RC': [] }
    for akey in tstr:
        tstr[akey] += '{:}'.format(dfa)

    wd = {}
    for aset in refd:
        wd[aset] = {}
        for arx in refd[aset]:
            wd[aset][arx] = {}

    dsplit = dfa.split('@')
    tdfa = dsplit[0]
    if len(dsplit) > 1:
        at_dfa = dsplit[1]

    tdir = scf_dir+'{:}/{:}_BH76/'.format(tdfa,dfa)

    tfl = tdir + 'BH76_total.yaml'
    td_scf_default = BH76_analysis(cdir=tdir)

    for aref in refs:

        #if dfa == aref:
        #    continue
        if (aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat):
            continue

        if aref == 'SCAN-FLOSIC':

            if '@' in dfa:
                tdir_scf = nrlmol_basis_dir+'{:}/{:}_BH76/'.format(tdfa, dfa)
                td_scf = BH76_analysis(cdir=tdir_scf)
            else:
                tdir_scf = './FLOSIC/{:}/{:}_BH76/'.format(tdfa,dfa)
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
                if '@' in dfa:
                    tdir = refs[aref]+'{:}/{:}@{:}_BH76/'.format(tdfa,tdfa,tref)
                else:
                    if dfa == aref:
                        tdir =  refs[aref]+'{:}/{:}_BH76/'.format(dfa,dfa)
                    else:
                        tdir = refs[aref]+'{:}/{:}@{:}_BH76/'.format(dfa,dfa,tref)

            tfl = tdir + 'BH76_total.yaml'
            BH76_analysis(cdir=tdir)
            td_nscf = yaml.load(open(tfl,'r'), Loader=yaml.Loader)

        for idict in range(2):

            if idict == 0:
                didx = 'RX'
                idx = 'BH76'
            elif idict == 1:
                didx = 'RC'
                idx = 'BH76RC'

            for asys in td_scf[didx]:
                tfe = (td_nscf[didx][asys]['Energy'] - refd[idx][asys])
                tde = (td_scf[didx][asys]['Energy'] - td_nscf[didx][asys]['Energy'])

                wd[idx][asys][aref] = {'FE': tfe, 'DE': tde}

    outd = {}
    for aset in wd:
        outd[aset] = {}
        mfe = 0.
        mafe = 0.
        mde = 0.
        made = 0.
        ucrt_fe = 0.
        ucrt_de = 0.
        npass = 0
        for arx in wd[aset]:
            FEs = [ wd[aset][arx][aref]['FE'] for aref in wd[aset][arx] ]
            DEs = [ wd[aset][arx][aref]['DE'] for aref in wd[aset][arx] ]
            FE_ucrt = max(FEs) - min(FEs)
            DE_ucrt = max(DEs) - min(DEs)
            if FE_ucrt < cutoff and DE_ucrt < cutoff:
                outd[aset][arx] = {
                    'FE': sum(FEs)/(1.*len(FEs)),
                    'FE ucrt': FE_ucrt,
                    'DE': sum(DEs)/(1.*len(DEs)),
                    'DE ucrt': DE_ucrt
                }
                npass += 1
            else:
                failed[aset].append(arx)
                fail_d[aset][arx-1] += 1
        FEs = np.array([ outd[aset][arx]['FE'] for arx in outd[aset] ])
        FE_ucrts = np.array([ outd[aset][arx]['FE ucrt'] for arx in outd[aset] ])
        DEs = np.array([ outd[aset][arx]['DE'] for arx in outd[aset] ])
        DE_ucrts = np.array([ outd[aset][arx]['DE ucrt'] for arx in outd[aset] ])

        fnpass = 1.*npass
        outd[aset]['Stats'] = {
            'MAFE': np.sum(np.abs(FEs))/fnpass,
            'MFE': np.sum(FEs)/fnpass,
            'FE ucrt': (np.sum(FE_ucrts**2)/fnpass)**(0.5),
            'FE ucrt savg': np.sum(FE_ucrts)/fnpass,
            'MADE': np.sum(np.abs(DEs))/fnpass,
            'MDE': np.sum(DEs)/fnpass,
            'DE ucrt': (np.sum(DE_ucrts**2)/fnpass)**(0.5),
            'DE ucrt savg': np.sum(DE_ucrts)/fnpass,
            'Npass': npass
        }
        for akey in outd[aset]['Stats']:
            if akey != 'Npass':
                outd[aset]['Stats'][akey] = float(outd[aset]['Stats'][akey])

        tex_str[aset] += dfa
        for astat in ['MFE', 'MDE', 'FE ucrt', 'DE ucrt', 'FE ucrt savg', 'DE ucrt savg']:
            tex_str[aset] += ' & {:.2f}'.format(outd[aset]['Stats'][astat])

        lchar = ''
        if idfa == len(dfas) - 1:
            lchar = '\\hline'
        tex_str[aset] += ' & {:} \\\\ {:} \n'.\
            format(outd[aset]['Stats']['Npass'],lchar)

    outd['Failed cutoff'] = failed.copy()

    ofl = out_dir + dfa + '_dc_dft_cut.yaml'
    with open(ofl,'w+') as tfl:
        yaml.dump(outd,tfl,Dumper=yaml.Dumper)

with open(out_dir + 'dc_dft_cut_summ.tex','w+') as tfl:
    for aset in tex_str:
        tfl.write(str_rep_rules(tex_str[aset]))

for aset in refd:
    ncommon = 0
    for arx in refd[aset]:
        if fail_d[aset][arx-1] == len(dfas):
            ncommon += 1
            print('{:} reaction {:} failed for all DFAs'.format(aset,arx))
        elif fail_d[aset][arx-1] == len(dfas)-1:
            print('{:} reaction {:} failed for all but one DFA'.format(aset,arx))

    if ncommon == 0:
        print('No {:} reactions failed for all DFAs'.format(aset))
