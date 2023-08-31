# -*- coding: utf-8 -*-

import  numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
# example 5.9
cas = 'annulus'
fluid = 'water'
patm = 101.325*1000.0
mp = 1.00
Cp = 4180
H = 250
SDRi = 11
SDRo = 13.5
kpi = 0.4
kpo = 1.5
kfill = 1.7
doi,doo = sdr_pipe(4,SDRo)
rb = doo/2
dii,dio = sdr_pipe(2,SDRi)
Ac_ann = pi*(doi**2 - dio**2)/4
Ac_inn = pi*dii**2/4
dh = doi - dio
fluid = 'Water'
Tmk = 300
rho = PropsSI('D','T',Tmk,'P',patm,fluid)
mu = PropsSI('V','T',Tmk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Tmk,'P',patm,fluid)
Prw = PropsSI('Prandtl','T',Tmk,'P',patm,fluid)
k = PropsSI('conductivity','T',Tmk,'P',patm,fluid)
ua = mp/(rho*Ac_ann)
Rea = rho*ua*dh/mu
uc = mp/(rho*Ac_inn)
Rec = rho*uc*dii/mu

hi = 2000
hoi = 2000
hoo = 2000
Rfi = 1/(pi*dii*hi)
Rfo = 1/(pi*dio*hoi)
Rfe = 1/(pi*doi*hoo)
Rpi = np.log(dio/dii)/(2*pi*kpi)
Rpo = np.log(doo/doi)/(2*pi*kpo)
Rfill = 0
R12 = Rfi + Rpi + Rfo
R1 = Rfe + Rpo + Rfill
Rb = R1
Ra = 4*R1*R12/(4*R1+R12)
eta = H/(mp*Cp*np.sqrt(Rb*Ra))
eta2 = H/(mp*Cp*np.sqrt(Rb*R12))
Rbs = Rb*eta/np.tanh(eta)
Rbs2 = Rb*(1 + eta**2/3)
Rbs3= Rb*eta2/np.tanh(eta2)
Rbs4 = Rb*(1 + eta2**2/3)
CCf = mp*Cp
Rs = 0

zB,TdB,TuB,z_B,TyB,Q = Calcul_TBeir_coaxial(R12,R1,CCf,H,Rs,cas)
thoB = TuB[0]
xsi = np.sqrt(Ra/(4*R1))
gam  = H/(CCf*R1)
eta = gam/(2*xsi)
zz = zB/H
if cas  == 'annulus':
    Thd2,Thu2 = Tcoax_annulus(zz,xsi,gam,eta)
else:
    Thd2,Thu2 = Tcoax_inner(zz,xsi,gam,eta)
R1t = R1*CCf/H
R12t = R12*CCf/H
if cas  == 'annulus':
    Thd3,Thu3,tho3,thi3 = Tcoax_annulus_flux(zz,R1t,R12t)
else:
    Thd3,Thu3,tho3,thi3 = Tcoax_inner_flux(zz,R1t,R12t)
TdB = Thd3/Thd3[0]
TuB = Thu3/Thd3[0]
#Rbsnn = R1*gam*(1+tho)/(1-tho)
p1 = plot(Thd2,zz,color = 'k',linestyle = '-',marker = 'x',markersize=11,label = 'Uniform Temperature')
p2 = plot(Thu2,zz,color = 'k',linestyle = '-',marker = 'x',markersize=11,label = '')
p3 = plot(TdB,zz,color = 'k',linestyle = '-',marker = 'o',markersize=9,label = 'Uniform Flux')
p4 = plot(TuB,zz,color = 'k',linestyle = '-',marker = 'o',markersize=9,label = '')
legend(fontsize = 14)
ax = gca()
ax.invert_yaxis()
grid(True,which='both')
#title('Coubes de pertes de charge  ')
fts = 16
ftst = 14
xlabel(' $\Theta$',fontname='Times new Roman',fontsize = fts)
ylabel('z',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
ecri = 0
eta = H/(mp*Cp*np.sqrt(Ra*Rb))
Rbs1 = Rb*eta/np.tanh(eta)
Rbsc= Rb*(1 + Ra*eta**2/(3*R12))
tho1 =Thu2[0]
tho2 =TuB[0]
Rbs2 = Rb*gam/2*(1 + tho1)/(1 - tho1)
print ('Rb UT= ',Rbs1,Rbs2)
Rbs3 = Rb*(1 + 1/(3*R1t*R12t))
Rbs4 = Rb*gam/2*(1 + tho2)/(1 - tho2)
print(Rbs3,Rbs4)
show()
