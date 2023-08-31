#coding: latin-1
# Exemple 4.4
from geothermal_md import *
import numpy as np
#
# data
#
aljr = 0.09 # m2/jr
alhr = aljr/24.0
rb = 0.08
ks = 2.0
To = 6
q1 = -1500.0
q2 = -3000.0
q3 = -7500.0
Rb = 0.052
t1 = 700
t2 = 796
tf = 800
Tf = 15
Fof = alhr*tf/rb**2
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
G1 = G_function_ils(Fof)
G2 = G_function_ils(Fof-Fo1)
G3 = G_function_ils(Fof-Fo2)
SCT = (q1*G1 + (q2-q1)*G2+ (q3-q2)*G3)/ks
L = (SCT + q3*Rb)/(To-Tf)
print ('L = ',L)
