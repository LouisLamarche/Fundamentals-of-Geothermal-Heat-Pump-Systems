from  numpy import *
from geothermal_md import *
# example 2.1
def Nu_ann_turb(Re,a,Pr):
    k1 = 1.07 + 900/ Re -0.63/(1+10*Pr)
    Res =  Re*((1+a**2)*log(a) + (1 - a**2))/((1-a)**2*log(a))
    f = (0.7817*log(Res)-1.5)**(-2)
    F1 = 0.75*a**(-0.17)
    F2 = 0.9 - 0.15*a**0.6
    Fa = (F1+F2)/(1+a)
    Nu = f/8* Re*Pr/(k1 + 12.7*sqrt(f/8)*(Pr**(2/3) - 1))*Fa
    return Nu
def Nu_ann_lam(a):
    Nu = 3.66 + (4 - 0.102/(a+0.02))*a**(0.04)
    return Nu
def Nu_circ_turb(Re,Pr):
    fi = (0.79*log(Re)-1.64)**(-2)
    Nu = fi/8*(Re-1000)*Pr/(1 + 12.7*sqrt(fi/8)*(Pr**(2/3) - 1))
    return Nu
def Nu_circ_lam():
    return 3.66

vp  = 0.4/1000   # m3/s
rho = 1020
mu = 0.0025
Pr = 20
k = 0.5
mp = vp*rho
#
kpi = 0.4
kpo = 0.8
doo = 0.04
to = 0.003
doi = doo - 2*to
dio = 0.025
ti = 0.002
dii = dio - 2*ti
#
# annular region
Ac_ann = pi*(doi**2 - dio**2)/4

Ac_inn = pi*dii**2/4
dh = doi - dio
a = dio/doi
u_ann = mp/(rho*Ac_ann)
Re_ann = rho*u_ann*dh/mu
if Re_ann < 2300:
    Nu_ann =  Nu_ann_lam(a)
elif Re_ann > 4000:
    Nu_ann =  Nu_ann_turb(Re_ann,a,Pr)
else:
    gam = (Re_ann - 2300)/(4000-2300)
    Nu_l =  Nu_ann_lam(a)
    Nu_t =  Nu_ann_turb(Re_ann,a,Pr)
    Nu_ann = (1 - gam)*Nu_l + gam*Nu_t
h_ann = Nu_ann*k/dh
#
# inner region
#
u_inn = mp/(rho*Ac_inn)
Re_inn = rho*u_inn*dii/mu
if Re_inn < 2300:
    Nu_inn =  Nu_circ_lam(a)
elif Re_inn > 4000:
    Nu_inn =  Nu_circ_turb(Re_ann,Pr)
else:
    gam = (Re_inn - 2300)/(4000-2300)
    Nu_l =  Nu_circ_lam(a)
    Nu_t =  Nu_circ_turb(Re_ann,Pr)
    Nu_inn = (1 - gam)*Nu_l + gam*Nu_t
h_inn = Nu_inn*k/dh
Rfi = 1/(pi*dii*h_inn)
Rfo = 1/(pi*dio*h_ann)
Rfe = 1/(pi*doi*h_ann)
Rpi = log(dio/dii)/(2*pi*kpi)
Rpo = log(doo/doi)/(2*pi*kpo)
R12 = Rfi + Rpi + Rfo
R1 = Rfe + Rpo
print('Re annulus = ',Re_ann)
print('h annulus = ',h_ann)
print('Re inner = ',Re_inn)
print('h inner = ',h_inn)
print('R12 = ',R12)
print('R1= ',R1)
q = -10000
H = 150
qp = q/H
Cp = 3800
CCf = mp*Cp
R1t = R1*CCf/H
R12t = R12*CCf/H
Rbs = R1*(1 + 1/(3*R1t*R12t))
rb = doo/2
t = 48
alj = 0.1
ks = 2.75
Fo = alj*2/rb**2
Rs = G_function(Fo)/ks
To = 8
Tb = To - qp*Rs
Tf = Tb - qp*Rbs
Tfo = Tf + q/(2*CCf)
dz = 0.05
zz = arange(0,1+dz,dz)
cas = 'annulus'
if cas  == 'annulus':
    Thd3,Thu3,tho3,thi3 = Tcoax_annulus_flux(zz,R1t,R12t)
else:
    Thd3,Thu3,tho3,thi3 = Tcoax_inner_flux(zz,R1t,R12t)
Tref = q/CCf
Tfob =  Tb-tho3*Tref
Td3 = Thd3*Tref + Tb
Tu3 = Thu3*Tref + Tb
print(Tfob,Tfo)
