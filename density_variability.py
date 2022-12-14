import yaml

from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion

conversion() # ensures FLOSIC errors are recomputed

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

dfas = ['LSDA','PBE','SCAN']
hf_dir = './results_aug-cc-pvqz/'
nrlmol_dir = './results_NRLMOL_cart/'
refs = {
    'LCwPBE': './results_aug-cc-pvqz/',
    'S50X': './scan_hybs/S50X/',
    'SCAN-FLOSIC': './FLOSIC/'
}

tstr = {'BH76': '', 'BH76RC': ''}
for aset in tstr:
    tstr[aset] += aset
    for aref in refs:
        tstr[aset] += ' & \\multicolumn{2}{c}{'+'{:}'.format(aref) +'}'
    tstr[aset] += ' \\\\ \nDFA'
    for aref in refs:
        tstr[aset] += ' & Mean & Mean Abs.'
    tstr[aset] += ' \\\\ \hline \n'

for idfa, dfa in enumerate(dfas):

    at_hf_default = BH76_analysis(cdir=hf_dir+'{:}/{:}@HF_BH76/'\
        .format(dfa,dfa),wrc=True)

    for aset in tstr:
        tstr[aset] += '{:}'.format(dfa)

    for aref in refs:

        devs = {'BH76': {'MDV': 0., 'MADV': 0.}, 'BH76RC': {'MDV': 0., 'MADV': 0.}}
        if aref == 'SCAN-FLOSIC':
            at_hf = BH76_analysis(cdir=nrlmol_dir+'{:}/{:}@HF_BH76/'\
                .format(dfa,dfa),wrc=True)
        else:
            at_hf = at_hf_default.copy()

        if aref in ['S25X','S50X']:
            tdir = refs[aref]+'{:}@{:}_BH76/'.format(dfa,aref)
        else:
            tdir = refs[aref]+'{:}/{:}@{:}_BH76/'.format(dfa,dfa,aref)

        if aref == 'SCAN-FLOSIC':
            at_ref = yaml.load(open(tdir + 'BH76_total.yaml','r'), Loader=yaml.Loader)
        else:
            at_ref = BH76_analysis(cdir=tdir)

        for aset in ['BH76','BH76RC']:

            if aset == 'BH76':
                idx = 'RX'
            elif aset == 'BH76RC':
                idx = 'RC'

            for arx in at_hf[idx]:
                tval = at_hf[idx][arx]['Energy'] - at_ref[idx][arx]['Energy']
                devs[aset]['MDV'] += tval
                devs[aset]['MADV'] += abs(tval)

            for akey in ['MDV','MADV']:
                devs[aset][akey] /= (1.*len(at_hf[idx].keys()))

            tstr[aset] += ' & {:.2f} & {:.2f}'\
                .format(devs[aset]['MDV'],devs[aset]['MADV'])

    lchar = ''
    if idfa == len(dfas)-1:
        lchar = '\\hline'
    for aset in tstr:
        tstr[aset] += ' \\\\ {:} \n'.format(lchar)

with open('./dens_var.tex','w+') as tfl:
    for aset in tstr:
        tfl.write(str_rep_rules(tstr[aset]))
