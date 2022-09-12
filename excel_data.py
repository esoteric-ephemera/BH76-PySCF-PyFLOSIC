import pandas as pd
from analysis import BH76_analysis

tfl = './BH76-DFA@FLOdata.xlsx'
td = pd.read_excel(tfl,usecols=[1,2,3,4,5],skiprows=1)

for dfa in td.columns[1:]:
    tmpd = {}
    for isys,asys in enumerate(td.iloc[:,0]):
        tmpd[asys] = float(td[dfa].iloc[isys])
    BH76_analysis(cdir='./',edict=tmpd,fprefix=dfa+'_',wrc=True)
