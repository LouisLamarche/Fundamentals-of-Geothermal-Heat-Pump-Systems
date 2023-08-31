
#
# Example 5.7written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md  import *
from  matplotlib.pyplot import *
#
cas =1   # 2 flow rate divided by 2
fluid = 'Water'
patm = 101.325*1000.0
als = 0.1       # m2/day
alhr = als/24.0 # m2/hr
rb = 0.15/2.0
ks = 2.0
To = 10.0
qb = -40.0
H = 100.0
#
# borehole
di,do = sdr_pipe(1.25,11)     #  SDR-11 1.25 po nominal
ro = do/2.0
ri = di/2.0
xc = rb/3
Vp = 0.8                # débit en l/s
Trefk = 30.0+273
rhof = PropsSI('D', 'T', Trefk, 'P', patm, fluid)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', patm, fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)

mp = Vp*rhof/1000.0
kp = 0.4       # conductivité du plastique
#
# Calcul de la résistance de la conduite
#
kg = 1.7   # coulis    mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
# pipe resistance
rcond = np.log(ro/ri)/(2*pi*kp)
Re = 4*mp/(pi*di*mu)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hf = (Nud*kf)/(di)
rconv = 1/(pi*di*hf)
Rp = rcond + rconv
sigma = (kg-ks)/(kg+ks)
if cas ==2:
    mp = mp/10
Rg,Ra = Rb_linesource(kg,ks,rb,ro,xc)
Rb = Rg + Rp/2
Ra = Ra + 2*Rp

rbb = rb/H
zob = 0.04
q1 = 1000.0
q2 = 1800.0
q3 = -200.0
q1p = q1/H
q2p = q2/H
q3p = q3/H
# intermediate time in hours
t1 = 1.0
t2 = 3.0
t3 = 5.0
dt = 0.1
t = np.arange(0,t3+dt,dt)
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
Fo3 = alhr*t3/rb**2
R1a = G_function(Fo3)
R2a = G_function(Fo3-Fo1)
R3a = G_function(Fo3-Fo2)
DT1 = (q1p*R1a+(q2p-q1p)*R2a+(q3p-q2p)*R3a)/ks
Tb1 = To - DT1  # Tb with CLS
print ('After 5 hours')
print ('Tb a) = ',Tb1)
Tf = Tb1 - q3p*Rb
q = q3p*H
CCf = mp*Cp
gam = H/CCf/Rb
Tfo = Tf + q/(2*CCf)
Tfi = Tf - q/(2*CCf)
print ('Tfo lineaire a) = ',Tfo)
print ('Tfi lineaire a) = ',Tfi)
#
x = gam/2.0
Rbe = Rb*x/np.tanh(x)
Tfe = Tb1 - q3p*Rbe
Tfoe = Tfe + q/(2*CCf)
Tfie = Tfe - q/(2*CCf)
print ('Tfo exp a) = ',Tfoe)
print ('Tfi exp a) = ',Tfie)
#

xsi = np.sqrt(Ra/(4*Rb))
eta =  gam/(2*xsi)
Rbeff = Rb*eta/np.tanh(eta)
Tfz = Tb1 - q3p*Rbeff
Tfoz = Tfz + q/(2*CCf)
Tfiz = Tfz - q/(2*CCf)
print ('Tfo Zeng a) = ',Tfoz)
print ('Tfi Zeng a) = ',Tfiz)
theol = (1-gam/2)/(1+gam/2)
theoe = np.exp(-gam)
theoz = (np.cosh(eta) - xsi*np.sinh(eta)) /(np.cosh(eta) + xsi*np.sinh(eta))

Tfol2 = Tb1 -q3p*Rb*gam*theol/(1-theol)
Tfoe2 = Tb1 -q3p*Rb*gam*theoe/(1-theoe)
Tfoz2 = Tb1 -q3p*Rb*gam*theoz/(1-theoz)
print('Tfo linear = ',Tfol2)
print('Tfo exp = ',Tfoe2)
print('Tfo Zeng = ',Tfoz2)


#
# Temperature at tjhe end of pulse 2

#
R1a = G_function(Fo2)
R2a = G_function(Fo2-Fo1)
DT1 = (q1p*R1a+(q2p-q1p)*R2a)/ks
Tb1 = To - DT1  # Tb with CLS
print ('After second pulse')
print ('Tb a) = ',Tb1)
Tf = Tb1 - q2p*Rb
q = q2p*H
CCf = mp*Cp
gam = H/CCf/Rb
Tfo = Tf + q/(2*CCf)
Tfi = Tf - q/(2*CCf)
print ('Tfo lineaire b) = ',Tfo)
print ('Tfi lineaire b) = ',Tfi)
#
x = gam/2.0
Rbe = Rb*x/np.tanh(x)
Tfe = Tb1 - q2p*Rbe
Tfoe = Tfe + q/(2*CCf)
Tfie = Tfe - q/(2*CCf)
print ('Tfo exp b) = ',Tfoe)
print ('Tfi exp b) = ',Tfie)
#

xsi = np.sqrt(Ra/(4*Rb))
eta =  gam/(2*xsi)
Rbeff = Rb*eta/np.tanh(eta)
Tfz = Tb1 - q2p*Rbeff
Tfoz = Tfz + q/(2*CCf)
Tfiz = Tfz - q/(2*CCf)
print ('Tfo Zeng b) = ',Tfoz)
print ('Tfi Zeng b) = ',Tfiz)
theol = (1-gam/2)/(1+gam/2)
theoe = np.exp(-gam)
theoz = (np.cosh(eta) - xsi*np.sinh(eta)) /(np.cosh(eta) + xsi*np.sinh(eta))

Tfol2 = Tb1 -q2p*Rb*gam*theol/(1-theol)
Tfoe2 = Tb1 -q2p*Rb*gam*theoe/(1-theoe)
Tfoz2 = Tb1 -q2p*Rb*gam*theoz/(1-theoz)
print('Tfo linear b= ',Tfol2)
print('Tfo exp b= ',Tfoe2)
print('Tfo Zeng b= ',Tfoz2)
flag_plot = False
if flag_plot:
#
#
    #
    #
    nh = len(t)
    qp = np.zeros(nh)
    Tb1v = np.zeros(nh)
    #
    # Calcululation of Tb
    #
    for i in range(0,nh):
        Fo = alhr*t[i]/rb**2
        if t[i] < t1:
            DT1 = q1p*G_function(Fo)/ks             # ICS
            qp[i] = q1p
        elif t[i] < t2:
            DT1 = (q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1))/ks
            qp[i] = q2p
        else:
            DT1 = (q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1)+(q3p-q2p)*G_function(Fo-Fo2))/ks
            qp[i] = q3p
        Tb1v[i] = To - DT1  # Tb with ICS
    print ('Tb a) = ',Tb1v[nh-1] )
    Tfv = Tb1v - qp*Rb
    q = qp*H
    Tfov = Tfv + q/(2*mp*Cp)
    Tfiv = Tfv - q/(2*mp*Cp)
    #print ('Tf = ',Tf[nh-1] )
    print ('Tfo = ',Tfov[nh-1] )
    #print ('Tfi = ',Tfi[nh-1] )
    Tfev = Tb1v - qp*Rbe
    q = qp*H
    Tfoev = Tfev + q/(2*mp*Cp)
    Tfiev = Tfev - q/(2*mp*Cp)
    #print ('Tf e= ',Tfe[nh-1] )
    print ('Tfo e= ',Tfoev[nh-1] )
    #print ('Tfi e= ',Tfie[nh-1] )
    Tfzv = Tb1v - qp*Rbeff
    q = qp*H
    Tfozv = Tfzv + q/(2*mp*Cp)
    Tfizv = Tfzv - q/(2*mp*Cp)
    #print ('Tf z = ',Tfz[nh-1] )
    print ('Tfo z = ',Tfozv[nh-1] )
    #print ('Tfi z= ',Tfiz[nh-1] )
    figure(1)
    plot(t,Tb1v,t,Tfov,t,Tfoev,t,Tfozv)
    legend(('Tb','Tfo lin','Tfo  exp','Tfo zeng'))
    figure(2)
    plot(t,Tb1v,t,Tfov,t,Tfiv)
    legend(('Tb','Tfo lin','Tfi  lin'))
    show()
