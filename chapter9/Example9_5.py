#coding: utf-8
#
# example 9.5 (Parameter estimation method)
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from scipy.optimize import curve_fit
import pandas as pd
#
# donnees du probleme
#
gam = 0.5772157
M = np.loadtxt("..\\data\\TRT_test1.txt")
t = M[:,0]      # time in hours
Tf = M[:,1]     # température Farenheight
Tc = (Tf-32.0)/1.8# température Celsius
nt = len(Tc)
qp = 55.4      # heat transfer per metre
nt = len(Tc)    # total number of points
To = 10.4  # Undisturbed temperature
CC = 2.5e6      # Volumetric heat capacity
CCv = np.array([0.9*2.5e6,2.5e6,1.1*2.5e6])
rb = 0.075      # rayon du puits en metres
rbv = [0.9*rb,rb,1.1*rb]
ksv = np.zeros((3,3))
Rbv = np.zeros((3,3))
diffk = np.zeros((3,3))
diffRb = np.zeros((3,3))
#
#
ks = 2.25        # assumed intial value
#
# first approach
#
Rb_nom = 0.1806
ks_nom = 2.2648
n_crit = 140       # From example 7_2
y = Tc[n_crit:nt]
ti = t[n_crit:nt]
ok = False
compt_max = 10
compt = 0
G_vect = np.vectorize(G_function_ils)

def T_theo(x, knew,Rbnew):
    return  To + qp*(x/knew + Rbnew)
#
#
Rb = Rb_nom
ks = ks_nom
for ic in range(0,3):
    CC = CCv[ic]
    for ir in range(0,3):
        rb = rbv[ir]
        ok = False
        compt = 0
        po = [ks,Rb] #  initial guess
        while not ok:
            als = ks/CC
            alhr = als*3600
            Fo = alhr*ti/rb**2
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
            ksv[ic,ir] = ks
            Rbv[ic,ir] = Rb

flag_plot = False
if flag_plot:
    ksnom = ksv[1,1]
    Rbnom = Rbv[1,1]
    CC = CCv[1]
    rb = rbv[1]
    als = ksnom/CC
    alhr = als*3600
    Fo = alhr*ti/rb**2
    x = G_vect(Fo)
    T_thr = T_theo(x,ksnom,Rbnom)
    p1 = plot(t,Tc,ti,T_thr,'*')
    #axis('equal')
    grid(True,which='both')
    fts = 16
    ftst = 14
    ylabel(r'$T_f\: (\circ$C)',fontsize = fts)
    xlabel('time(min)',fontsize = fts)
    xticks(fontsize=ftst)
    yticks(fontsize=ftst)
    show()

diffk = abs(ksv - ks_nom)/ks_nom
diffRb = abs(Rbv - Rb_nom)/Rb_nom
indd  = ['r$_b$ = 6.75 cm','r$_b$ = 7.5 cm','r$_b$ = 8.25 cm']
coll = ['(rho C_p) = 2.25 MJ/m3-K','(rho C_p) = 2.5 MJ/m3-K','(rho C_p) = 2.75 MJ/m3-K']
def f3(x):
    return '%.3f' % x
def f2(x):
    return '%.3f' % x
def f1(x):
    return '%.1f' % x

df1 = pd.DataFrame(ksv,columns = coll,index = indd)
print(df1.to_string(index=True,formatters=[f2,f2,f2]).replace('nan',''))
df2 = pd.DataFrame(diffk*100,columns = coll,index = indd)
print(df2.to_string(index=True,formatters=[f1,f1,f1]).replace('nan',''))
df1 = pd.DataFrame(Rbv,columns = coll,index = indd)
print(df1.to_string(index=True,formatters=[f2,f2,f2]).replace('nan',''))
df2 = pd.DataFrame(diffRb*100,columns = coll,index = indd)
print(df2.to_string(index=True,formatters=[f1,f1,f1]).replace('nan',''))
