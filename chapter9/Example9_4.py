#coding: utf-8
#
# example de la methode de la ligne source pour évaluer la conductivité du sol
#
import numpy as np
from geothermal_md import *
from conversion_md import *
from  matplotlib.pyplot import *
import pandas as pd
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
DTo = 0.5
CC = 2.5e6      # Volumetric heat capacity
CCv = np.array([0.9*2.5e6,2.5e6,1.1*2.5e6])
rb = 0.075      # rayon du puits en metres
rbv = [0.9*rb,rb,1.1*rb]
ntests = 9;
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
cas = 'b'
if cas == 'a':
    G_vect = np.vectorize(G_function_ils)
else:
    G_vect = np.vectorize(G_function)
for ic in range(0,3):
    CC = CCv[ic]
    for ir in range(0,3):
        rb = rbv[ir]
        ok = False
        compt = 0
        while not ok:
            als = ks/CC
            alhr = als*3600
            Fo = alhr*ti/rb**2
            x = G_vect(Fo)
            p = np.polyfit(x,y,1)
            m = p[0]
            b = p[1]
            ksn = qp/m
            if (abs(ksn-ks)/ks) < 0.001:
                ok = True
            else:
                ks = ksn
                compt = compt + 1
                if compt > compt_max:
                    err = 1
                    ok = True
            Rb = (b-To)/qp
            ksv[ic,ir] = ksn
            Rbv[ic,ir] = Rb
diffk = abs(ksv - ks_nom)/ks_nom
diffRb = abs(Rbv - Rb_nom)/Rb_nom
indd  = ['rb = 6.75 cm','rb = 7.5 cm','rb = 8.25 cm']
coll = ['(rho C_p) = 2.25 MJ/m3-K','(rho C_p) = 2.5 MJ/m3-K','(rho C_p) = 2.75 MJ/m3-K']
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