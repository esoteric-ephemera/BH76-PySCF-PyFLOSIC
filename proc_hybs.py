from glob import glob
import numpy as np
import yaml

from analysis import BH76_analysis

tlist = glob('./*_EXX_BH76')
slname = tlist[0].split('/')[-1].split('_')[0]
tlist = [slname + '_BH76'] + tlist

outs = np.zeros((len(tlist),3))

for idir,adir in enumerate(tlist):
    if idir == 0:
        outs[idir,0] = 0.0
    else:
        outs[idir,0] = float(adir.split('/')[-1].split('_')[1])/100.
    BH76_analysis(cdir=adir + '/',wrc=True)
    td = yaml.load(open(adir + '/BH76_total.yaml','r'),Loader=yaml.Loader)
    outs[idir,1] = td['Stats']['MAD']
    outs[idir,2] = td['RC Stats']['MAD']

outs = outs[np.argsort(outs[:,0])]

with open('hyb_mads.csv','w+') as tfl:
    tfl.write('x, BH76 MAD (kcal/mol), BH76RC MAD (kcal/mol)\n')
    for idir in range(len(tlist)):
        tfl.write('{:}, {:}, {:}\n'.format(*outs[idir]))
