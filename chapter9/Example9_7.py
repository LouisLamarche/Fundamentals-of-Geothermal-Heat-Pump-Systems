#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
#
#
# Example 9_6
#
#
gamm = np.euler_gamma
qp = 66
rb = 0.126/2.0
To = 12.0
CC = 1.92e6
ksa = 2.5   # assumed
al = ksa/CC*3600.0  # m2/hr
alm = ksa/CC*60.0  # m2/hr
xmax = np.log(1 + al*50/(5*rb**2))
print ('xmax = ',xmax)
p1 = [3.5,18.5]
p2 = [0.5,13]
m1 =(p1[1]-p2[1])/(p1[0]-p2[0])
k1 = qp/(4*pi*m1)
print ('ks = ',k1)
#
# b)
#
xmaxb = np.log(5*rb**2/alm)
print ('xmaxb = ',xmaxb)
p1 = [6.0,21.5]
p2 = [8.,25.4]
m2 = (p1[1]-p2[1])/(p1[0]-p2[0])
k2 = qp/(4*pi*m2)
print ('ks2 = ',k2)
b = p1[1] - m2*p1[0]
als = k1/CC
alm = als*60
Rb = (b-To)/qp - (np.log(4*alm/rb**2)-gamm)/(4*pi*k1)
print ('Rb = ',Rb)

To = 12.0
CC = 1.92e6
rb = 0.126/2.0
qp = 66       # flux de chaleur par metre
t1 = 3000.0  # temps ou on arrete l'injection de chaleur
#
M = np.loadtxt("..\\data\\TRT_test2.txt")
t = M[:,0]     # temps en minutes
Tfi = M[:,1]
Tfo = M[:,2]
Tf2 = (Tfi+Tfo)/2.0
Tf = M[:,3]
nt = len(t)
n_inj = np.sum(t < t1)     # nombre de points d'injection
n_res = nt - n_inj       # nombre de points recovery
Tf_inj = Tf[0:n_inj]
Tf_res = Tf[n_inj:nt]
t_inj = t[0:n_inj]
t_res = t[n_inj:nt]

n1 = 10        # Minimum number of points to take out
nf = 10        # Minimum number of points to keep for the regressoion
ntot = n_res - nf - n1 + 1
#
ksolr = np.zeros(ntot)
for i in range(0,ntot):
    j1 = n1+i
    tr = t_res[j1:n_res]
    y = Tf_res[j1:n_res]
    ok = False
    compt = 0
    x = np.log(tr/(tr-t1))
    p = np.polyfit(x,y,1)
    m = p[0]
    ksolr[i] = qp/(4.0*pi*m)
nlast = 10
ksr = np.mean(ksolr[ntot-nlast:ntot-1])
print (ksr)
xk = n1+ np.arange(0,ntot)
nnn = 10
plot(xk[0:ntot:nnn],ksolr[0:ntot:nnn],'ko')

gam = 0.5772157
n1 = 10        # Minimum number of points to take out
nf = 10         # Minimum number of points to keep for the regressoion
ntot = n_inj - nf - n1 + 1
#
ksoli = np.zeros(ntot)
Rbi = np.zeros(ntot)
for i in range(0,ntot):
    j1 = n1+i
    tr = t_inj[j1:n_res]
    y = Tf_inj[j1:n_res]
    x = np.log(tr)
    p = np.polyfit(x,y,1)
    m = p[0]
    b = p[1]
    ksoli[i] = qp/(4.0*pi*m)
    almin =   ksoli[i]/CC*60.0 # m2/min
    Rbi[i] = (b-To)/qp - (np.log(4*almin/rb**2)-gam)/(4*pi*ksoli[i])

ksi = np.mean(ksoli[ntot-nlast:ntot-1])
Rbi = np.mean(Rbi[ntot-nlast:ntot-1])
print (ksi,Rbi)
xk = n1+ np.arange(0,ntot)
plot(xk[0:ntot:nnn],ksoli[0:ntot:nnn],'k+')
legend(('Restitution','Injection'),fontsize = 16)
xticks(fontsize= 14)
yticks(fontsize=14)
show()


