# -*- coding: utf-8 -*-

import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from CoolProp.CoolProp import *

patm = 101.325*1000.0
data = np.loadtxt('..\\data\\mmc1b.txt')
t1 = data[:,0]*60 # in seconds
Tin = data[:,1]
Tout = data[:,2]
Tf = data[:,3]
ne = len(t1)
nx = 30
xx = t1[0:ne:nx]
yy = Tf[1:ne:nx]
nn = len(xx)
rb = 0.126/2.0
L = 18.3
do = 0.03340
di = 0.02733
ro = do/2.0
ri = di/2
xc = 0.053/2
ks = 2.82
kg = 0.73
mp = 0.197
q = 1056
qp = q/L
CCg = 3.84e6
CCs = 1.92e6
als = ks/CCs
alg = kg/CCg
hinv = 1000
kp = 0.4
Trefk = 305
fluid = 'Water'
Cpf = PropsSI('Cpmass','T',Trefk,'P',patm,fluid)
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
CCf = rhof*Cpf

# pipe resistance
Re = 4*mp/(pi*di*mu)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hf = (Nud*kf)/(di)
Rp = 1/(4*pi*kp)*np.log(do/di) + 1/(2*pi*di*hinv)
J = 10
z = np.array([xc,-xc])  # pipes coordinates
Rbmp,Ramp = Rb_multipole(kg,ks,rb,ro,Rp*2,J,z)
Rbc1,Rac = Rb_linesource(kg,ks,rb,ro,xc)

Rgrout = 0.105
flag_increase = False
cas = 'a'
if cas == 'a':
    req = rb/np.exp(Rgrout*2*pi*kg)
    algrout = alg*(rb**2-req**2)/(rb**2-2*ro**2)
    gamma = np.sqrt(als/algrout)
    kt = kg/ks
    CCfeq = CCf*2*(ri/req)**2
    ret = req/rb
    nu = CCfeq*ret**2/(2*CCs)
    if flag_increase:
        nu = 0.1365
    Rt = 2*pi*ks*Rp
    fic = 'fig13_1a.png'
elif cas == 'b':
    req = np.sqrt(2)*ro
    kgrout_eq = np.log(rb/req)/(2*pi*Rgrout)
    algrout = kgrout_eq/CCg
    gamma = np.sqrt(als/algrout)
    kt =  kgrout_eq/ks
    CCfeq = CCf*2*(ri/req)**2
    ret = req/rb
    nu = CCfeq*ret**2/(2*CCs)
    if flag_increase:
        nu = 0.164
    Rt = 2*pi*ks*Rp
    fic = 'fig13_1b.png'
else:
    fic = 'fig13_1c.png'
if cas == 'a' or cas == 'b':
    Tn1 = np.zeros(nn)
    Tn2 = np.zeros(nn)
    Tn3 = np.zeros(nn)


    for i in range(0,nn):
        Fo = als*xx[i]/rb**2
        gg1 = G_function_st(Fo,ret,gamma,kt,nu,Rt)
        Tn1[i] = 22 + qp*gg1/ks
        gg2 = G_function(Fo)
        Tn2[i] = 22+ qp*(gg2/ks + Rgrout + Rp)
        gg3 = G_function_ils(Fo)
        Tn3[i] = 22+ qp*(gg3/ks + Rgrout + Rp)

    x = np.log(xx)
    p1 = plot(x,yy,color = 'k',linewidth = 2,label = 'Sand-Box data')
    p2 = plot(x,Tn1,color = 'k',linestyle = 'None',marker = 'o',markersize=8,label = 'short-time response')
    p3 = plot(x,Tn2,color = 'k',linestyle = '-',marker = 'x',markersize=8,label = 'ICS')
    p4 = plot(x,Tn3,color = 'k',linestyle = '-',marker = '+',markersize=8,label = 'ILS')
    legend(fontsize = 14)
    ax = gca()
    #dx = ax.get_xlim()
    #dy = ax.get_ylim()
    #ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
    grid(True,which='both')
    fts = 16
    ftst = 14

    xlabel(' log(t)',fontname='Times new Roman',fontsize = fts)
    ylabel(' $T_f (\circ$ C)',fontsize = fts,fontname='Times new Roman')


    xticks(fontsize=ftst)
    yticks(fontsize=ftst)
    tight_layout()
else:
    Rb1 = 0.105
    CCf = rhof*Cpf
    req1 = rb/np.exp(Rb1*2*pi*kg)
    ale1 = alg*(rb**2-req1**2)/(rb**2-2*ro**2)
    CCe1 = CCg*(rb**2-2*ro**2)/(rb**2-req1**2)
    gam1 = np.sqrt(als/ale1)
    kt1 = kg/ks
    ret1 = req1/rb
    nu1 = 0.1353
    CCfeq1 = nu1*CCs*(rb/ri)**2
    MF1 =  CCfeq1/CCf
    R1 = 2*pi*ks*Rp
    req2 = np.sqrt(2)*ro
    kg2 = np.log(rb/req2)/(2*pi*Rb1)
    alg2 = kg2/CCg
    gam2 = np.sqrt(als/alg2)
    kt2 = kg2/ks
    ret2 = req2/rb
    nu2 = 0.1643
    CCfeq2 = nu2*CCs*(rb/ri)**2
    MF2 =  CCfeq2/CCf
    Tna = np.zeros(nn)
    Tnb = np.zeros(nn)

    for i in range(0,nn):
        Fo = als*xx[i]/rb**2
        gg = G_function_st(Fo,ret1,gam1,kt1,nu1,R1)
        Tna[i] = 22+ qp*gg/ks
        gg = G_function_st(Fo,ret2,gam2,kt2,nu2,R1)
        Tnb[i] = 22+ qp*gg/ks
    x = np.log(xx)
    p1 = plot(x,yy,color = 'k',linewidth = 2,label = 'Sand-Box data')
    p2 = plot(x,Tna,color = 'k',linestyle = 'None',marker = 'o',markersize=8,label = '$r_{eq}$ ( Eq 13.34)')
    p3 = plot(x,Tnb,color = 'k',linestyle = 'None',marker = 'x',markersize=8,label = '$r_{eq}$ ( Eq 13.30)')
    le = legend(fontsize = 14)

    ax = gca()
    #dx = ax.get_xlim()
    #dy = ax.get_ylim()
    #ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
    grid(True,which='both')
    #title('Coubes de pertes de charge  ')
    fts = 16
    ftst = 14
#    xlabel(' log(t)',fontname='Times new Roman',fontsize = fts)
#    ylabel(' Tf($\circ$ C)',fontsize = fts,fontname='Times new Roman')
    xlabel(' log(t)',fontname='Times new Roman',fontsize = fts)
    ylabel(' $T_f (\circ$ C)',fontsize = fts,fontname='Times new Roman')

    xticks(fontsize=ftst)
    yticks(fontsize=ftst)
    tight_layout()
show()
