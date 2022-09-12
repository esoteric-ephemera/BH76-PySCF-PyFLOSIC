import yaml
from analysis import BH76_analysis,rdir
from os import path

ostr = ',|--,MD,--,--|,|--,MAD,--,--|,|--,RMSD,--,--|\n'
ostr += 'DFA'
for i in range(3):
    ostr += ', DFT, DFT@HF, DFT@HF(TS), DFT@LCwPBE'
ostr += '\n'
dfa_l = ['HF','LSDA','PBE','SCAN','r2SCAN','LCwPBE']

for dfa in dfa_l:

    dft_fl = dfa + '_BH76/BH76_total.yaml'
    if not path.isfile(dft_fl):
        BH76_analysis(cdir=dfa+'_BH76/')
    dft_d = yaml.load(open(dft_fl,'r'),Loader=yaml.Loader)

    if dfa == 'HF':
        ostr += '{:}, {:}, , , , {:}, , , , {:}, , , \n'.format(dfa,\
            dft_d['Stats']['MD'],dft_d['Stats']['MAD'],dft_d['Stats']['RMSD'])
        continue

    hfdft_fl = dfa + '@HF_BH76/BH76_total.yaml'
    if not path.isfile(hfdft_fl):
        BH76_analysis(cdir=dfa+'@HF_BH76/')

    hfdft_d = yaml.load(open(hfdft_fl,'r'),Loader=yaml.Loader)

    refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
        Loader=yaml.Loader)

    en_d = {}

    for indx in refs:
        td = refs[indx]['Stoich']
        for aspec in td:
            if td[aspec] > 0.0:
                en_d[aspec] = hfdft_d['SP'][aspec]
            else:
                en_d[aspec] = dft_d['SP'][aspec]

    BH76_analysis(cdir=dfa+'_BH76/',edict=en_d,fprefix=dfa+'_HF_TS_')

    tnm = dfa+'_BH76/'+dfa+'_HF_TS_BH76_total.yaml'
    hf_ts_d = yaml.load(open(tnm,'r'),Loader=yaml.Loader)

    if dfa == 'LCwPBE':
        lc_dft_d = {'Stats': {'MD': '', 'MAD': '', 'RMSD': ''}}
    else:
        lc_dft_fl = dfa + '@LCwPBE_BH76/BH76_total.yaml'
        if not path.isfile(lc_dft_fl):
            BH76_analysis(cdir=dfa+'@LCwPBE_BH76/')

        lc_dft_d = yaml.load(open(lc_dft_fl,'r'),Loader=yaml.Loader)

    ostr += dfa
    for stat in ['MD','MAD','RMSD']:
        for adict in [dft_d,hfdft_d,hf_ts_d,lc_dft_d]:
            ostr += ', {:}'.format(adict['Stats'][stat])
    ostr += '\n'

with open('BH76_comp.csv','w+') as tfl:
    tfl.write(ostr)
