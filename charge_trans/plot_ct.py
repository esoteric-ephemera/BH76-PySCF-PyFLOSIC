import numpy as np
import matplotlib.pyplot as plt

colors = ['darkblue','darkorange','tab:green']
labels = [['FCl...CH$_3$',.5,-.04],['F...H$_2$',.63,.02],
['Cl...H$_2$',.1,-.016]]
tdat = np.genfromtxt('./tsct.csv',delimiter=',',skip_header=1,\
    skip_footer=1)

fig, ax = plt.subplots(figsize=(6,4))

tmpmat = np.vstack([tdat[:,0], np.ones(tdat.shape[0])]).T
ixl = np.linspace(0.,1.,2000)

for icol in range(tdat.shape[1]-1):
    ax.scatter(tdat[:,0],tdat[:,icol+1],color=colors[icol],\
        label=labels[icol][0])
    slp, ict = np.linalg.lstsq(tmpmat,tdat[:,icol+1],rcond=None)[0]
    print(labels[icol],-ict/slp)
    ax.plot(ixl,slp*ixl+ict,color=colors[icol],linestyle=':')
    ax.annotate(labels[icol][0],(labels[icol][1],labels[icol][2]),\
        color=colors[icol],fontsize=12)

ax.set_xlim(0.,1.)
ax.set_xlabel('$x$',fontsize=12)
ax.set_ylabel('$\\Delta N$',fontsize=12)

ax.hlines(0.,0.,1.,color='k',linestyle='-',linewidth=1)

#ax.legend(fontsize=12)
#plt.show()
plt.savefig('./charge_trans_ts.pdf',dpi=600,bbox_inches='tight')
