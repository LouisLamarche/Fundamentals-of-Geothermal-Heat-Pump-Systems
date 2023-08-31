# -*- coding: utf-8 -*-

import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from CoolProp.CoolProp import *

patm = 101.325*1000.0
rb = 0.065
do = 0.03340
di = 0.02733
ro = do/2.0
ri = di/2
xc = 0.053/2
ks = 2.5
kg = 1.25
CCg = 2.e6
CCs = 2.2e6
CCs2 = 2.5e6/2
als = ks/CCs
alhr = als*3600
alg = kg/CCg
Trefk = 305
fluid = 'Water'
Cpf = PropsSI('Cpmass','T',Trefk,'P',patm,fluid)
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
CCf = rhof*Cpf
mp = 0.197
kp = 0.4
# pipe resistance
Re = 4*mp/(pi*di*mu)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hinv = (Nud*kf)/(di)
Rp = 1/(4*pi*kp)*np.log(do/di) + 1/(2*pi*di*hinv)
J = 10
z = np.array([xc,-xc])  # pipes coordinates
Rbmp,Ramp = Rb_multipole(kg,ks,rb,ro,Rp*2,J,z)
Rbc1,Rac = Rb_linesource(kg,ks,rb,ro,xc)
Rgrout = Rbmp
Rgrout2 = 0.93*Rbmp
req = rb/np.exp(Rgrout*2*pi*kg)
req2 = rb/np.exp(Rgrout2*2*pi*kg)
algrout = alg*(rb**2-req**2)/(rb**2-2*ro**2)
algrout2 = alg*(rb**2-req2**2)/(rb**2-2*ro**2)
gamma = np.sqrt(als/algrout)
kt = kg/ks
CCf = CCs
CCfeq = CCf*2*(ri/req)**2
ret = req/rb
nu = CCfeq*ret**2/(2*CCs)
Rp = 0
Rt = 2*pi*ks*Rp
tt = np.arange(0,100,1/10)
nn = len(tt)
To = 11
t1 = 160
t2 = t1 + 8
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
q1p = -55
q2p = 0
q3p = -60
Tn1 = np.zeros(nn)
Tn2 = np.zeros(nn)
qp = np.zeros(nn)
for i in range(0,nn):
    Fo = alhr*tt[i]/rb**2
    if tt[i] < t1:
        DT1 = q1p*G_function_st(Fo,ret,gamma,kt,nu,Rt)/ks             # ICS
        DT2 = q1p*G_function_ils(Fo)/ks             # ICS
        qp[i] = q1p
    elif tt[i] < t2:
        DT1 = (q1p*G_function_st(Fo,ret,gamma,kt,nu,Rt)/ks+(q2p-q1p)*G_function_st(Fo-Fo1,ret,gamma,kt,nu,Rt))/ks
        qp[i] = q2p
    else:
        DT1 = (q1p*G_function_st(Fo,ret,gamma,kt,nu,Rt)+(q2p-q1p)*G_function_st(Fo-Fo1,ret,gamma,kt,nu,Rt)\
            +(q3p-q2p)*G_function_st(Fo-Fo2,ret,gamma,kt,nu,Rt))/ks
        qp[i] = q3p
    Tn1[i] = To - DT1  # Tb with ICS
    Tn2[i] = To - DT2  # Tb with ICS

x = np.log(tt)
A = np.zeros((nn,3))
A[:,0] = tt*60
A[:,1] = Tn1
A[:,2] = Tn2
np.savetxt('test_trt2.txt',A)
p2 = plot(x,Tn1,color = 'k',linestyle = 'None',marker = 'o',markersize=8,label = 'short-time response')
p3 = plot(x,Tn2,color = 'k',linestyle = '-',marker = 'None',markersize=8,label = 'short-time response')
legend(fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14
xlabel(' log(t)',fontname='Times new Roman',fontsize = fts)
ylabel(' $T_f (\circ$ C)',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()
