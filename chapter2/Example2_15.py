import numpy as np
from utilities_md import *
# exemple 2.15
data = np.loadtxt("..\\data\\Temp_Mtl.txt")  # loads in  kWatts
Tamb1 = data[:,1]    # external temperatures
nt = len(Tamb1)
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])
Ktot =  31.0
gains = 225.0
Ti = 22.0
Te = Ti - gains/Ktot
ib = 4
Tmax = max(Tamb1)
Tmin = min(Tamb1)
DT = Tmax-Tmin
Tm = (Tmax+Tmin)/2
nbin = int(np.ceil(DT/ib))
Tminn = np.floor(Tmin)
Tmaxn = Tminn+nbin*ib
#
# solution using Histogram function
#
Nb,Tedge = np.histogram(Tamb1,bins = nbin ,range = (Tminn,Tmaxn))
Nba = np.zeros(nbin)
DT = np.zeros(nbin)
Nbb = np.zeros((nbin,2))
n1a = 0
n2a = 0
n1b = 0
n2b = 0
for i in range(0,8760):
    T = Tamb1[i]
    ii = cherche_index2(T,Tedge)
    ih = i % 24
    if ih > 8 & ih < 17:
        Nbb[ii,1] = Nbb[ii,1] + 1
    else:
        Nbb[ii,0] = Nbb[ii,0] + 1
file = open('exe2_8.txt','w')
Qbina = np.zeros(nbin)
Qbinb = np.zeros(nbin)
for i in range(0,nbin):
    strr = '%d' % Tedge[i] +  ':'+ '%d' % Tedge[i+1]  + '&' + '%d' % Nbb[i,0] + '&' + '%d' % Nbb[i,1] +   '\\\ ' + '\n'
    file.write(strr)
Ntot1  = sum(Nbb[:,0])
Ntot2  = sum(Nbb[:,1])
strr = 'Total  '+   '&' + '%d' % Ntot1+ '&' +'%d' % Ntot2 + '\\\ ' + '\n'
file.write(strr)
file.close()
file = open('exe2_8b.txt','w')
for i in range(0,nbin):
    Tbin = Tedge[i] + ib/2.0
    DT[i] = max(Te - Tbin,0)
    Qbina[i] = Ktot*Nbb[i,0]*max(Te - Tbin,0)/1000.0
    Qbinb[i] = Ktot*Nbb[i,1]*max(Te - Tbin,0)/1000.0
    strr = '%d' % Tedge[i] +  ':'+ '%d' % Tedge[i+1]  + '&' + '%.1f' % DT[i] + '&' + '%.1f' % Qbina[i] +  '&' +'%.1f' %  Qbinb[i] + '\\\ ' + '\n'
    file.write(strr)
Qtota  = sum(Qbina)
Qtotb  = sum(Qbinb)
strr = 'Total  '+   '&'+   '&' + '%.1f' % Qtota+ '&' +'%.1f' % Qtotb + '\\\ ' + '\n'
file.write(strr)
file.close()
Q = Qtota +Qtotb
print ('c) Heat losses using BIN method are  ')
print ('%.2f' % Q + ' kWh')
