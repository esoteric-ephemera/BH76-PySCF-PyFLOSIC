from itertools import product
from os import sys

v1 = range(0,4,1)
v2 = range(0,8,1)
v3 = range(0,5,1)

def vfail(v1,v2,v3):
    print("class FLO: This code requires pyscf >= 1.5.2; your version {:}.{:}.{:}".format(v1,v2,v3))
    #sys.exit()
    return

def vpass(v1,v2,v3):
    print("Version passed; your version {:}.{:}.{:}".format(v1,v2,v3))
    return

for pv in product(v1,v2,v3):

    vmaj,vmin,vfix = pv

    if vmaj > 1:
        vpass(vmaj,vmin,vfix)
    elif vmaj == 1:
        if vmin > 5:
            vpass(vmaj,vmin,vfix)
        elif vmin == 5 and vfix >= 2:
            vpass(vmaj,vmin,vfix)
        else:
            vfail(vmaj,vmin,vfix)
    else:
        vfail(vmaj,vmin,vfix)
