#coding: utf-8
#
# exemple 11.6
#
import  numpy as np
from geothermal_mod import *
from conversion_mod import *
from CoolProp.CoolProp import *


cas = 'b'
rhoc = 1000
muc = 0.0014
Prc = 10
Cpc = 4200
kc = 0.57
#
rhoh = 1020
muh = 0.0025
Prh = 20
Cph = 3900
kh = 0.5
do = 0.10
di = 0.06
mpc = 1
mph = 0.5
Cc = mpc*Cpc
Ch = mph*Cph
L = 10
Tci = 7.0
Thi = 40.0


def Nu_ann_turb(Re,a,Pr):
    k1 = 1.07 + 900/ Re -0.63/(1+10*Pr)
    Res =  Re*((1+a**2)*log(a) + (1 - a**2))/((1-a)**2*log(a))
    f = (0.7817*log(Res)-1.5)**(-2)
    Fa = 0.75*a**(-0.17)
    Nu = f/8* Re*Pr/(k1 + 12.7*sqrt(f/8)*(Pr**(2/3) - 1))*Fa
    return Nu
def Nu_ann_lam(a):
    Nu = 3.66 + 1.2*a**(-0.8)
    return Nu
def Nu_circ_turb(Re,Pr):
    fi = (0.79*log(Re)-1.64)**(-2)
    Nu = fi/8*(Re-1000)*Pr/(1 + 12.7*sqrt(fi/8)*(Pr**(2/3) - 1))
    return Nu
def Nu_circ_lam():
    return 3.66


def Calcul_NTU(eff,Cr):
    if Cr < 1:
        NTU = 1/(Cr-1)*np.log((eff-1)/(eff*Cr-1))
    else:
        NTU = eff/(1.0-eff)
    return NTU

def Calcul_eff(NTU,Cr):

    if Cr < 1:
        eff = (1-np.exp(-NTU*(1-Cr)))/(1-Cr*np.exp(-NTU*(1-Cr)))
    else:
        eff = NTU/(1+NTU)
    return eff

#
#
if cas == 'a':
    Acc = pi*(do**2 - di**2)/4
    Ach = pi*di**2/4
    dc = do - di
    dh = di
    a = di/do
    #
    uc = mpc/(rhoc*Acc)
    Rec = rhoc*uc*dc/muc
    if Rec < 2300:
        Nuc =  Nu_ann_lam(a)
    elif Rec > 4000:
        Nuc =  Nu_ann_turb(Rec,a,Prc)
    else:
        gam = (Rec - 2300)/(4000-2300)
        Nu_l =  Nu_ann_lam(a)
        Nu_t =  Nu_ann_turb(Rec,a,Prc)
        Nuc = (1 - gam)*Nu_l + gam*Nu_t
    hc = Nuc*kc/dc
    #
    uh = mph/(rhoh*Ach)
    Reh = rhoh*uh*dh/muh
    if Reh < 2300:
        Nuh =  Nu_circ_lam()
    elif Reh > 4000:
        Nuh =  Nu_circ_turb(Reh,Prh)
    else:
        gam = (Reh - 2300)/(4000-2300)
        Nu_l =  Nu_circ_lam()
        Nu_t =  Nu_circ_turb(Reh,Prh)
        Nuh = (1 - gam)*Nu_l + gam*Nu_t
    hh = Nuh*kh/dh
else:
    Ach = pi*(do**2 - di**2)/4
    Acc = pi*di**2/4
    dh = do - di
    dc = di
    a = di/do
    #
    uc = mpc/(rhoc*Acc)
    Rec = rhoc*uc*dc/muc
    if Rec < 2300:
        Nuc =  Nu_circ_lam(a)
    elif Rec > 4000:
        Nuc =  Nu_circ_turb(Rec,Prc)
    else:
        gam = (Rec - 2300)/(4000-2300)
        Nu_l =  Nu_circ_lam()
        Nu_t =  Nu_circ_turb(Rec,Prc)
        Nuc = (1 - gam)*Nu_l + gam*Nu_t
    hc = Nuc*kc/dc
    #
    uh = mph/(rhoh*Ach)
    Reh = rhoh*uh*dh/muh
    if Reh < 2300:
        Nuh =  Nu_ann_lam(a)
    elif Reh > 4000:
        Nuh =  Nu_ann_turb(Reh,a,Prh)
    else:
        gam = (Reh - 2300)/(4000-2300)
        Nu_l =  Nu_ann_lam(a)
        Nu_t =  Nu_ann_turb(Reh,a,Prh)
        Nuh = (1 - gam)*Nu_l + gam*Nu_t
    hh = Nuh*kh/dh
print(Rec,Reh)
print(hh,hc)
# inner region
#
Ah = pi*di*L
Ac = Ah
Rh = 1/(hh*Ah)
Rc = 1/(hc*Ac)
Rt = Rh + Rc
U = hh*hc/(hh+hc)
UA = 1/Rt
UAb = U*Ah
Cmin = min(Cc,Ch)
Cmax = max(Cc,Ch)
qmax = Cmin*(Thi - Tci)
Cr = Cmin/Cmax
NTU = UA/Cmin
eff = Calcul_eff(NTU,Cr)
print(eff)
q = eff*qmax
print(q)
Tco = Tci + q/Cc
Tho = Thi - q/Ch
print(Tco,Tho)

