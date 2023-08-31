import numpy as np
from utilities_md import *
from matplotlib.pylab import *
from tabulate import tabulate
import pandas as pd
data = np.loadtxt("..\data\Temperature.txt",skiprows =1 )  # loads in  kWatts
Tamb1 = data[:,1]    # external temperatures
#Tamb1 = 32 + Tamb1
nt = len(Tamb1)
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])
COP = 3
Ktot =  1000.0
gainsA = 3000
gainsB = 6000
TiB = 22.0
TiA = 18.0
TeA = TiA - gainsA/Ktot
TeB = TiB - gainsB/Ktot
ib = 3
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
DTA = np.zeros(nbin)
DTB = np.zeros(nbin)
Nbb = np.zeros((nbin,6),dtype=int)
n1a = 0
n2a = 0
n1b = 0
n2b = 0
for i in range(0,8760):
    T = Tamb1[i]
    ii = cherche_index2(T,Tedge)
    ih = i % 24
    if ih > 8 and ih < 17:
        Nbb[ii,2] = Nbb[ii,2] + 1
        if T < TeB:
            n1b = n1b + 1
        else:
            n2b = n2b + 1
    else:
        Nbb[ii,0] = Nbb[ii,0] + 1
        if T < TeA:
            n1a = n1a + 1
        else:
            n2a = n2a + 1
totala = np.sum(Nbb,axis=1)
Nb = Nbb[:,2]*5/7
Na  = Nbb[:,2] - Nb.astype(int)
totalb =  Nbb[:,0] + Na + Nb.astype(int)
Nbb[:,1] = Nbb[:,0] + Na
Nbb[:,3] = Nb.astype(int)
Nbb[:,4] = totala
Nbb[:,5] = totalb
coll = ['A','An','B','Bn','total','total']
indd = []
Tmi = int(Tmin)
for i in range(0,nbin):
    Tma = Tmi + ib - 1
    pt = '%d' % Tmi+':'+'%d' % Tma
    indd.append(pt)
    Tmi = Tma + 1
df1 = pd.DataFrame(Nbb,columns = coll,index = indd)
print(df1)
Qbina = np.zeros(nbin)
Qbinb = np.zeros(nbin)
Pbina = np.zeros(nbin)
Pbinb = np.zeros(nbin)
for i in range(0,nbin):
    Tbin = Tedge[i] + ib/2.0
    Qbina[i] = Ktot*Nbb[i,0]*max(TeA - Tbin,0)/1000.0
    Qbinb[i] = Ktot*Nbb[i,1]*max(TeB - Tbin,0)/1000.0
    if Tbin <= -10:
        Pbina[i] = Qbina[i]
        Pbinb[i] = Qbinb[i]
    else:
        Pbina[i] = Qbina[i]/COP
        Pbinb[i] = Qbinb[i]/COP
Qtota  = sum(Qbina)
Qtotb  = sum(Qbinb)
Q = Qtota +Qtotb
Ptota  = sum(Pbina)
Ptotb  = sum(Pbinb)
P = Ptota +Ptotb
print ('A: Q = ' + '%.2f' % Qtota + ' kWh, P = '+ '%.2f' % Ptota + ' kWh')
print ('B: Q = ' + '%.2f' % Qtotb + ' kWh, P = '+ '%.2f' % Ptotb + ' kWh')
print ('Q = ' + '%.2f' % Q + 'P = ' + '%.2f' % P + ' kWh')
df1['Qa'] = Qbina
df1['Qb'] = Qbinb
print(df1)

