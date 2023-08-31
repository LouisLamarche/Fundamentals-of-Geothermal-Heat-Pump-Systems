# -*- coding: utf-8 -*-

import  numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
# example 5.10
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
xsi = np.sqrt(Ra/(4*R1))
gam  = H/(CCf*R1)
eta = gam/(2*xsi)
zB,TdB,TuB,z_B,TyB,Q = Calcul_TBeir_coaxial(R12,R1,CCf,H,0,cas)
zz = zB/H
if cas  == 'annulus':
    Thd2,Thu2 = Tcoax_annulus(zz,xsi,gam,eta)
else:
    Thd2,Thu2 = Tcoax_inner(zz,xsi,gam,eta)
sca = 50
reo = rb
Lc = sca*reo
data1 = np.loadtxt('Tup_comsolai.txt')
n1 = len(data1)
ni = 0
To = 0
fact = 1
z1a= data1[:,0]
z2a = np.zeros(n1)
for i in range(0,n1):
    z2a[i] = z1a[n1 - i - 1]*Lc
zu_a = z1a[ni:n1-1]*Lc
Tu_a = data1[ni:n1-1,1]
Tu_comsol = fact*(data1[ni:n1-1,1]) + To
data2 = np.loadtxt('Tdown_comsolai.txt')
n2 = len(data2)
zd_a = data2[ni:n2-1,0]*Lc
Td_a = data2[ni:n2-1,1]
Td_comsol = fact*(data2[ni:n2-1,1]) + To


R1t = R1*CCf/H
R12t = R12*CCf/H
if cas  == 'annulus':
    Thd3,Thu3,tho3,thi3 = Tcoax_annulus_flux(zz,R1t,R12t)
else:
    Thd3,Thu3,tho3,thi3 = Tcoax_inner_flux(zz,R1t,R12t)
TdB3 = Thd3/Thd3[0]
TuB3 = Thu3/Thd3[0]
Fo = 100
ks = 2.5
Rs = G_function_ils(Fo)/ks
qp = -50
To = 0
Tb = To - qp*Rs
Tbi = Tb
ni = 0
zB,TdB,TuB,z_B,TyB,Q = Calcul_TBeir_coaxial(R12,R1,CCf,H,Rs,cas)
thoBeier = TuB[0]
tho_Tunif =Thu2[0]
tho_qunif =TuB3[0]
q = qp*H
Tfi = To - q/(CCf*(1-thoBeier))
Tfix = To + q/(CCf*(thoBeier-1 ))
Tfo = thoBeier*(Tfi-To) + To
print('Beirs Tin = ',Tfi)
print('Beirs Tout = ',Tfo)
Tf = TyB*(Tfi - To) + To
Td_B = TdB*(Tfi - To) + To
Tu_B = TuB*(Tfi - To) + To
DT  = q/(CCf*(1-tho_Tunif))
Td2 = -Thd2*DT + Tb
Tu2 = -Thu2*DT + Tb
print('Unif Tb, Tin = ',Td2[0])
print('Unif Tb, Tout = ',Tu2[0])

Tref = q/CCf
Td3 = -Thd3*Tref + Tb
Tu3 = -Thu3*Tref + Tb
DTq  = q/(CCf*(1-tho_qunif))
Td3b = -TdB3*DTq + Tb
Tu3b = -TuB3*DTq + Tb
print('Unif q, Tin = ',Td3[0])
print('Unif q, Tout = ',Tu3[0])
nfu = len(Tu_comsol) - 1
nfd = len(Td_comsol) - 1
nx = 40
Tunn = Tu_comsol[0:nfu:nx]
Tdnn = Td_comsol[0:nfd:nx]
zun = zu_a[0:nfu:nx]/250
zdn = zd_a[0:nfd:nx]/250
#Rbsnn = R1*gam*(1+tho)/(1-tho)
p1 = plot(Td2,zz,color = 'k',linestyle = '-',marker = 'x',markersize=11,label = 'Uniform Temperature ')
p2 = plot(Tu2,zz,color = 'k',linestyle = '-',marker = 'x',markersize=11,label = '')
p3 = plot(Td3,zz,color = 'k',linestyle = '-',marker = 'o',markersize=9,label = 'Uniform Flux')
p4 = plot(Tu3,zz,color = 'k',linestyle = '-',marker = 'o',markersize=9,label = '')
p5 = plot(Td_B,zz,color = 'k',linestyle = '-',marker = '+',markersize=9,label = 'Beier')
p6 = plot(Tu_B,zz,color = 'k',linestyle = '-',marker = '+',markersize=9,label = '')
p7 = plot(Tunn,zun,color = 'k',linestyle = ':',marker = 'v',markersize=9,label ='Comsol')
p8 = plot(Tdnn,zdn,color = 'k',linestyle = ':',marker = 'v',markersize=9,label = '')
legend(fontsize = 14)
ax = gca()
ax.invert_yaxis()
grid(True,which='both')
#title('Coubes de pertes de charge  ')
fts = 16
ftst = 14
xlabel(' T($^\circ$ C)',fontname='Times new Roman',fontsize = fts)
ylabel('z',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
eta = H/(mp*Cp*np.sqrt(Ra*Rb))
Rbs1 = Rb*eta/np.tanh(eta)
Rbsc= Rb*(1 + Ra*eta**2/(3*R12))
show()
