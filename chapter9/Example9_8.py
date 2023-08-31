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

To = 12.0
CC = 1.92e6
rb = 0.126/2.0
gam = 0.5772157
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
#
# Conductivity Estimation in the recovery
#
ks = 2.0        # assumed intial value
als = ks/CC
almin = als*60.0
Fomin  = 7
t_critr = t1 + Fomin*rb**2/almin     # tmin  hypothetique
n_critr = np.sum(t_res < t_critr)    # nombre de point minimum à enlever
y = Tf_res[n_critr:n_res]
tr = t_res[n_critr:n_res]
nte = len(y)
ok = False
compt_max = 10
compt = 0
x = np.zeros(nte)
cas = 'a'
if cas == 'a':
    G_vect = np.vectorize(G_function_ils)
else:
    G_vect = np.vectorize(G_function)

def T_theo_rec(x, knew):
    return  To + qp*x/knew

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
print ('k=',ks)
T_thr = T_theo_rec(x,ks)
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
flag_plot = True
if flag_plot:
    figure(1)
    p1 = plot(t_res,Tf_res,tr,T_thr,'*')
    p2 = plot(t_inj,Tf_inj,ti,T_thi,'*')
    #axis('equal')
    grid(True,which='both')
    fts = 16
    ftst = 14
    ylabel(r'$T_f\: (\circ$C)',fontsize = fts)
    xlabel('time(min)',fontsize = fts)
    xticks(fontsize=ftst)
    yticks(fontsize=ftst)

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
T_thin = T_theo(x,ks,Rb)

rcParams.update({'figure.autolayout': True})
flag_plot = True
if flag_plot:
    figure(2)
    p1 = plot(t_res,Tf_res,tr,T_thr,'*')
    p2 = plot(t_inj,Tf_inj,ti,T_thin,'*')
    #axis('equal')
    grid(True,which='both')
    fts = 16
    ftst = 14
    ylabel(r'$T_f\: (\circ$C)',fontsize = fts)
    xlabel('time(min)',fontsize = fts)
    xticks(fontsize=ftst)
    yticks(fontsize=ftst)
    show()



