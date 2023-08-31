import  numpy as np
from geothermal_md import *
from scipy.optimize import fsolve

CC = 2.0e6
ks = 2.5
als = ks/CC
alhr = als*3600
alj = alhr*24
rb = 0.08
d = 3
r1 = 0.62
r2 = d - r1
r3 = d
t = 48
qp = -50
DTb_mes = 6.9
DT1_mes = 1.0
def fct(par):
    knew,Cnew = par
    alnew = knew/Cnew
    alhrnew = alnew*3600
    Xb = rb/(2*np.sqrt(alhrnew*t))
    X1 = r1/(2*np.sqrt(alhrnew*t))
    X2 = r2/(2*np.sqrt(alhrnew*t))
    X3 = r3/(2*np.sqrt(alhrnew*t))
    Ib = I_function(Xb)
    I1 = I_function(X1)
    I2 = I_function(X2)
    I3 = I_function(X3)
    y1  = -qp*(Ib+I3)/(2*pi*knew) - DTb_mes
    y2 = -qp*(I1+I2)/(2*pi*knew) - DT1_mes
    return(y1,y2)
param = fsolve(fct,(2,2e6))
ksol = param[0]
alhrn = ksol/param[1]*3600
print ('k sol = ' + '%.1f' % ksol)
print ('CC = ' + '%.1f' % (alhrn*1e4))
# ils2
Fo = alhr*t/rb**2
Fo1 = alhr*t/r1**2
Fo2 = alhr*t/r2**2
Fo3 = alhr*t/r3**2
Gb = G_function_ils(Fo)
G1 = G_function_ils(Fo1)
G2 = G_function_ils(Fo2)
G3 = G_function_ils(Fo3)
DTb = -qp*(Gb+G3)/(ks)
DT1 = -qp*(G1+G2)/(ks)
print ('DT (rb) = ' + '%.4f' % DTb)
print ('DT (r1) = ' + '%.4f' % DT1)
