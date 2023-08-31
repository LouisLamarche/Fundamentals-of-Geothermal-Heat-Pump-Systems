#coding: utf-8
#
# Example 5.6 written by Louis Lamarche 22 sept 2017
#
from numpy import *
from geothermal_md import *
from  matplotlib.pyplot import *
Rp = 0.08
ro = 0.026/2.0
ri = 0.0228/2.0
kp = 0.42
Rp = log(ro/ri)/(2*pi*kp)
kg = 1.5   # grout
ks = 2.1   # soil
rhos = 2213
rb = 0.075
xc = 0.015
Qf = 1.0e-4
mp = 1000*Qf
Cp_f = 4200
Cp_s = 1100;
als = ks/(rhos*Cp_s)
alhr = als*3600
CCf = mp*Cp_f
H = 150
q = 7500
qp = q/H
To = 10
#
# line source
#
Fo =  alhr*100/rb**2
Tb = To - qp/ks*G_function_ils(Fo)
#
# line source
#
sigma = (kg-ks)/(kg+ks)
Rg,Rga = Rb_linesource(kg,ks,rb,ro,xc)
print(Rg,Rga)
xct = xc/rb
Rp  = 0.05
Rb1 = Rg + Rp/2.0
Ra1 = Rga + 2*Rp
Rb2 = 1/(4*pi*kg)*(log(rb**2/(2*xc*ro))-sigma*log(1 - xct**4)) + Rp/2
Ra2 = 1/(pi*kg)*(log(2*xc/ro)+ sigma*log((1 + xct**2)/(1 - xct**2)))  + 2*Rp
print(Rb1,Ra1,Rb2,Ra2)
Tf = Tb - qp*Rb1
Tfol = Tf + q/(2*CCf)
Tfil = Tf - q/(2*CCf)
#
# exponential
gam = H/CCf/Rb1
x = gam/2.0
Rbe = Rb1*x/tanh(x)
Tfe = Tb - qp*Rbe
Tfoe = Tfe + q/(2*CCf)
Tfie = Tfe - q/(2*CCf)
#
# Zeng's profile
#
xsi = sqrt(Ra1/(4*Rb1))
eta =     gam/(2*xsi)
Rbz = Rb1*eta/tanh(eta)
Tfz = Tb - qp*Rbz
Tfoz = Tfz + q/(2*CCf)
Tfiz = Tfz - q/(2*CCf)
print('Tfo linear = ',Tfol)
print('Tfo exp = ',Tfoe)
print('Tfo Zeng = ',Tfoz)
theol = (1-gam/2)/(1+gam/2)
theoe = exp(-gam)
theoz = (cosh(eta) - xsi*sinh(eta)) /(cosh(eta) + xsi*sinh(eta))
print(theol)
print(theoe)
print(theoz)
Tfol2 = Tb -qp*Rb1*gam*theol/(1-theol)
Tfoe2 = Tb -qp*Rb1*gam*theoe/(1-theoe)
Tfoz2 = Tb -qp*Rb1*gam*theoz/(1-theoz)
print('Tfo linear = ',Tfol2)
print('Tfo exp = ',Tfoe2)
print('Tfo Zeng = ',Tfoz2)

zx,Td,Tu,zy,Ty = Compute_Taxial(Ra1,Rb1,CCf,H)
eta = H/(CCf*sqrt(Ra1*Rb1))
xsi = sqrt(Ra1/Rb1)/2.0
Tdd = (cosh(eta*(1 - zx/H)) + xsi*sinh(eta*(1 - zx/H)))/(cosh(eta) + xsi*sinh(eta))
Tuu = (cosh(eta*(1 - zx/H)) - xsi*sinh(eta*(1 - zx/H)))/(cosh(eta) + xsi*sinh(eta))
Ran = 4*Rb1
zx,Td,Tu,zy,Ty2 = Compute_Taxial(Ran,Rb1,CCf,H)
#
# Compute real temperature
#
#
data = loadtxt('LWT_TConstantn.txt')
z = data[:,2]
T1a = data[:,100]
T1b = data[:,50]
data = loadtxt('EWT_TConstantn.txt')
z = data[:,2]
T2a = data[:,100]
T2b = data[:,50]
data = loadtxt('Tbn.txt')
ti = data[6,0]
t1 = data[6:223,0] - ti
Tbv = data[6:223,1]
nt = len(t1)
Tbn = Tbv[100]
#Tfi =
Th1a = T1a - Tbn
Th2a = T2a - Tbn
Th1b = T1b - Tbn
Th2b = T2b - Tbn
nz = len(z) - 1
The1a  = Th1a/Th1a[nz]
The2a  = Th2a/Th1a[nz]
The1b  = Th1b/Th1b[nz]
The2b  = Th2b/Th1b[nz]
zt = z/-150;
zt2 = zy/150;

T_thz = Ty*(Tfiz-Tb)+Tb
T_the = Ty2*(Tfie-Tb)+Tb
nx = len(zx)
nnx = 2*nx
Tyl = zeros(nnx)
for i in range(0,nx):
    Td = Tfil + (Tfol-Tfil)/(2*H)*zx[i]
    Tu = Tfil + (Tfol-Tfil)/(2*H)*(2*H-zx[i])
    Tyl[i] = Td
    Tyl[nnx-i - 1] = Tu
nnz = 2*nz
T_com = zeros(nnz)
z_com = zeros(nnz)
for i in range(0,nz):
    T_com[nz - i - 1] = T1a[i]
    T_com[nz + i] = T2a[i]
    z_com[nz - i - 1] = z[i]
    z_com[nz + i] = z[i]
Tfoc = T_com[nnz-1]
Tfic = T_com[0]
theoc =  (Tfoc -   Tb)/(Tfic - Tb)
print(theol)
print(theoe)
print(theoz)
print(theoc)
errl = abs(theol - theoc)/theoc
erre = abs(theoe - theoc)/theoc
errz = abs(theoz - theoc)/theoc
print(errl)
print(erre)
print(errz)
ftst  = 16
fts = 14
plot(T_com,z_com,T_thz,-zy,'x',T_the,-zy,'o',Tyl,-zy,'-+')
legend(('Comsol','Zeng profile','exponential profile','linear profile'),fontsize = fts,frameon=False)
xticks(fontsize=ftst)
yticks(fontsize=ftst)
show()