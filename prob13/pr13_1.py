#coding: latin-1
# Exemple 4.4
from geothermal_md import *
import numpy as np
#
# data
#
aljr = 0.09 # m2/jr
alhr = aljr/24.0
algjr = 0.045 # m2/jr
alghr = algjr/24.0
ks = 2.0
To = 6
q1 = -1500.0
q2 = -3000.0
q3 = -7500.0
rb = 0.08 # feet
do = 0.043 # feet
di = do
ro = do/2.0
ri = di/2
xc = 0.5*rb # feet
kg =  2
alsjr = 0.09
alhr = alsjr/24
als = alhr/3600
algjr = 0.052
alghr = algjr/24
alg = alghr/3600
Rb,Ra = Rb_linesource(kg,ks,rb,ro,xc)
Rg = Rb
req = rb/np.exp(Rg*2*pi*kg)
ale = alg*(rb**2-req**2)/(rb**2-2*ro**2)
CCge = kg/ale
CCs = ks/als
gam = np.sqrt(als/ale)
kt = kg/ks
CCf = 4.2e6
CCfe = CCf*2*(ri/req)**2
ret = req/rb
nu = CCfe*ret**2/(2*CCs)
R = 0
t1 = 700
t2 = 796
tf = 800
Tf = 15
Fof = alhr*tf/rb**2
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
G1 = G_function_ils(Fof)
G2 = G_function_ils(Fof-Fo1)
G3 = G_function_st(Fof-Fo2,ret,gam,kt,nu,R)
R1 = (G1 - G2)/ks
R2 = (G2  + ks*Rb - G3)/ks
R3 = G3/ks
SCT = q1*R1 + q2*R2+ q3*R3
L = SCT/(To-Tf)
print ('L = ',L)
