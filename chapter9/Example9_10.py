#coding: utf-8
#
# exemple  9.7
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from scipy.optimize import curve_fit
#
#

CC = 1.92e6
To = 8.22
#To = 7.894
CC = 2.25e6
rb = 0.057
gam = np.euler_gamma
qp = 66.42       # flux de chaleur par metre
t1 = 4861.0  # temps ou on arrete l'injection de chaleur
#
M = np.loadtxt("..\\data\\TRT_test3.txt")
ts = M[:,0]     # temps en secondes
t = ts/60.0
#t = ts/3600.0
Tfi = M[:,4]   # temune inérature en Celsius
Tfo = M[:,8]   # température en Celsius
Tf = (Tfi+Tfo)/2.0
DTin = Tfi - To
DTout = Tfo - To
p = -0.99999
cas = 'linear'
if cas == 'p-linear':
    Tf = To +  p*(abs(DTin)**(p+1) -abs(DTout)**(p+1))/((1+p)*(abs(DTin)**(p) -abs(DTout)**(p)))
nt = len(t)
n_inj = np.sum(t < t1)     # nombre de points d'injection
n_res = nt - n_inj       # nombre de points recovery
Tf_inj = Tf[0:n_inj]
Tf_res = Tf[n_inj:nt]
t_inj = t[0:n_inj]
t_res = t[n_inj:nt]
#
# Conductivity Estimation in the recovery
#
ks = 2.0        # assumed intial value
als = ks/CC
almin = als*60.0
Fomin  = 8
t_critr = t1 + Fomin*rb**2/almin     # tmin  hypothetique

n_critr = np.sum(t_res < t_critr)    # nombre de point minimum à enlever
y = Tf_res[n_critr:n_res]
tr = t_res[n_critr:n_res]
nte = len(y)
compt_max = 10
compt = 0
x = np.zeros(nte)
G_vect = np.vectorize(G_function_ils)

def T_theo_rec(x, knew):
    return  To + qp*x/knew

# second approach
n1 = 10        # Minimum number of points to take out
nf = 3000         # Minimum number of points to keep for the regressoion
ntot = nte - nf - n1 + 1
#
ksol = np.zeros(ntot)
ksolp = np.zeros(ntot)
nkept = np.zeros(ntot)
for i in range(0,ntot):
    j1 = n1+i
    tr = t_res[j1:n_res]
    y = Tf_res[j1:n_res]
    nkept[i] = len(y)
    ok = False
    compt = 0
    x = np.log(tr/(tr-t1))
    p = np.polyfit(x,y,1)
    m = p[0]
    ksolp[i] = qp/(4.0*pi*m)
    while not ok:
        als = ks/CC
        almin = als*60
        Fo1 = almin*t1/rb**2
        Fo = almin*tr/rb**2
        x = G_vect(Fo) - G_vect(Fo-Fo1)
        param,resn = curve_fit(T_theo_rec,x,y)
        ksn = param[0]
        if (abs(ksn-ks)/ks) < 0.001:
            ok = True
        else:
            ks = ksn
            compt = compt + 1
            if compt > compt_max:
                err = 1
                ok = True
    ksol[i] = ks
nlast = 20
ksf = np.mean(ksol[ntot-nlast:ntot-1])
print (ksf)
xk = n1+ np.arange(0,ntot)
figure(1)
nnn = 100
p1 = plot(nkept[0:ntot:nnn],ksol[0:ntot:nnn],'k+')
p2 = plot(nkept[0:ntot:nnn],ksolp[0:ntot:nnn],'ko')
legend(('Parameter estimation ','Slope method'),fontsize = 16)
xlim(max(nkept), min(nkept))
xticks(fontsize= 14)
yticks(fontsize=14)
ylabel(r'k ( W/m-K)',fontsize = 14)
xlabel('Number of points kept',fontsize = 14)
T_thr = T_theo_rec(x,ks)
figure(2)
#p1 = plot(t_res,Tf_res,tr,T_thr,'*')
#
# Estimation de Rb during 'injection
#
als = ks/CC
almin = als*60
t_criti = Fomin*rb**2/almin # tmin  hypothetique
n_criti = sum(t_inj < t_criti) # nombre de point minimum à enlever
y = Tf_inj[n_criti:n_inj]
ti = t_inj[n_criti:n_inj]
nte = len(y)
x = np.zeros(nte)
if cas == 'a':
    G_vect = np.vectorize(G_function_ils)
else:
    G_vect = np.vectorize(G_function)
Fo = almin*ti/rb**2
x = G_vect(Fo)
ok = False
compt_max = 10
compt = 0
Rb = 0.05       # valeur inititale
def T_theo_inj(x, Rbnew):
    return To + qp*(x/ks + Rbnew)
param,resn = curve_fit(T_theo_inj,x,y)
Rb = param[0]
print ('Rb = ',Rb)
T_thi = T_theo_inj(x,Rb)

#
# c) estimation of k and Rb during injection
#
def T_theo(x, knew,Rbnew):
    return  To + qp*(x/knew + Rbnew)
ok = False
compt = 0
po = [ks,Rb] #  initial guess
while not ok:
    als = ks/CC
    almin = als*60
    Fo = almin*ti/rb**2
    x = G_vect(Fo)
    params,resn = curve_fit(T_theo,x,y,po)
    ksn = params[0]
    Rbn = params[1]
    if (abs(ksn-ks)/ks) < 0.001 and (abs(Rbn-Rb)/Rb) < 0.001:
        ok = True
    else:
        ks = ksn
        Rb = Rbn
        compt = compt + 1
        if compt > compt_max:
            err = 1
            ok = True
print ('ks c) =  ',ks)
print ('Rb c) =  ',Rb)


rcParams.update({'figure.autolayout': True})
#cas = 2
#if cas ==1:
p1 = plot(t_res,Tf_res,tr,T_thr,'*')
#else:
p2 = plot(t_inj,Tf_inj,ti,T_thi,'*')
#axis('equal')
grid(True,which='both')
fts = 16
ftst = 14
ylabel(r'$T_f\: (\circ$C)',fontsize = ftst)
xlabel('time(min)',fontsize = fts)
xticks(fontsize=ftst)
yticks(fontsize=ftst)
show()



