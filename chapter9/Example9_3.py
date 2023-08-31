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
CCv = [0.9*2.5e6,2.5e6,1.1*2.5e6]
rb = 0.075      # rayon du puits en metres
rbv = [0.9*rb,rb,1.1*rb]
Tov = [To-DTo,To,To+DTo]
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
n_crit = 140       # From example 9_2
y = Tc[n_crit:nt]
ti = t[n_crit:nt]
x = np.log(ti)
for ic in range(0,3):
    CC = CCv[ic]
    for ir in range(0,3):
#        rb = rbv[ir]
        To = Tov[ir]
        p = np.polyfit(x,y,1)
        m1 = p[0]
        b1 = p[1]
        ks = qp/(4.0*pi*m1)
        alhr =   ks/CC*3600.0 # m2/hr
        Rb = (b1-To)/qp - (np.log(4*alhr/rb**2)-gam)/(4*pi*ks)
        ksv[ic,ir] = ks
        Rbv[ic,ir] = Rb
diffk = (ksv - ks_nom)/ks_nom
diffRb = abs(Rbv - Rb_nom)/Rb_nom

print(Rbv)
indd  = ['To = 9.9 C','To = 10.4 C','To = 10.9 C']
coll = ['(rho C_p) = 2.25 MJ/m3-K','(rho C_p) = 2.5 MJ/m3-K','(rho C_p) = 2.75 MJ/m3-K']
def f2(x):
    return '%.3f' % x
def f1(x):
    return '%.1f' % x

df = pd.DataFrame(Rbv,columns = coll,index = indd)
print(df.to_string(index=True,formatters=[f2,f2,f2]).replace('nan',''))
df2 = pd.DataFrame(diffk*100,columns = coll,index = indd)
df3 = pd.DataFrame(diffRb*100,columns = coll,index = indd)
print(df3.to_string(index=True,formatters=[f1,f1,f1]).replace('nan',''))
