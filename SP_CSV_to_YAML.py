from os import path,system
from analysis import BH76_analysis

def conversion():
    wdict = {}

    with open('./FLOSIC/BH76data_flosic.csv','r') as infl:
        for irow,arow in enumerate(infl):
            tmp = arow.strip().split(',')
            if irow == 0:
                for i in range(1,len(tmp)):
                    wdict[tmp[i].strip()] = {}
                dfas = list(wdict.keys())
                Ndfa = len(dfas)
                continue
            mol = tmp[0]
            for i in range(Ndfa):
                wdict[dfas[i]][mol] = float(tmp[1+i])

    for dfa in dfas:
        tdir = './FLOSIC/'+dfa+'_BH76/'
        if not path.isdir(tdir):
            system('mkdir -p '+tdir)
        BH76_analysis(cdir=tdir,edict=wdict[dfa],wrc=True)
    return

if __name__ == "__main__":

    conversion()
