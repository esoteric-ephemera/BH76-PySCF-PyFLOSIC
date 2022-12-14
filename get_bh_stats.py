import yaml

td = yaml.load(open( './BH76_ref_energies.yaml','r'),Loader=yaml.Loader)
flen = 1.*len(td.keys())

mbh = 0.
mabh = 0.
var = 0.
nthresh = 0
smax = -1e20
smin = 1e20
amax = 0.
amin = 1e20

for akey in td:
    mbh += td[akey]['Ref']
    mabh += abs(td[akey]['Ref'])
    var += td[akey]['Ref']**2
    smax = max(smax,td[akey]['Ref'])
    smin = min(smin,td[akey]['Ref'])
    amax = max(amax,abs(td[akey]['Ref']))
    amin = min(amin,abs(td[akey]['Ref']))

mbh /= flen
mabh /= flen
var /= flen
rmsd = var**(0.5)
stddev = max(var - mbh**2,0.)**(0.5)
print('Mean(BH) = {:.2f}'.format(mbh))
print('Mean(|BH|) = {:.2f}'.format(mabh))
print('RMSD(BH) = {:.2f}'.format(rmsd))
print('STDDEV(BH) = {:.2f}'.format(stddev))
print('MIN(BH) = {:.2f}'.format(smin))
print('MAX(BH) = {:.2f}'.format(smax))
print('MIN(|BH|) = {:.2f}'.format(amin))
print('MAX(|BH|) = {:.2f}'.format(amax))
