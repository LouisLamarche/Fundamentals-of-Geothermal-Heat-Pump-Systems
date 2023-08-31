#
# exemple de la methode de la ligne source pour évaluer la conductivité du sol
#
import numpy as np
from geothermal_md import *
from finance_md import *
from  matplotlib.pyplot import *
from scipy.optimize import newton

disc = 0.08                # taux d'actualisation
Nccv = 15                  # nombre d,annéees  pour le calcul de la VAN
inf_comb = 0.04         # inflation sur combustible
inf_ord = 0.025         # inflation sur combustible
Npret = 10
inte = 0.04
cout_base = 8000.0
mise_depart = 0.10
deduct_impot = 0.20
Investa = 45000.0       # cout total
cout_caria = 4000.0
vra = 2000
Investb = 38000.0       # cout total
cout_carib = 4800.0
vrb = 1000
maintenancea = np.zeros(Nccv)
maintenanceb = np.zeros(Nccv)
nn = 4
maintenanceb[nn-1:Nccv:nn] = 500.0


def calcul_npv(t_act,dpa,paiementi,econ_combi,vr,mi,maintenance):
    glycol_cost = 0
    for i in range(0,Nccv):
        cost = maintenance[i]*(1+inf_ord)**i/(1+t_act)**(i+1)
        glycol_cost = glycol_cost + cost
    vrp = present_value(vr,Nccv,0,t_act)
    VAN = -dpa  - paiementi*pwf(Npret,0,t_act) +  econ_combi*pwf(Nccv,inf_comb,t_act) \
        - glycol_cost + vrp + deduct_impot*mi*pwf_int(Npret,Nccv,inte,t_act)
    return VAN


def calcul_npvv(t_act,dpa,paiementi,econ_combi,vr,mi,maintenance):

    num = 0
    den = dpa
    paiement = paiementi
    j = 0
    econ_combustible = econ_combi
    interet = mi*inte
    cap_paye = paiement-interet
    balance  = mi - cap_paye
    costs = paiement
    gains = econ_combustible + deduct_impot*interet
    economies = gains - costs
    va = present_value(economies,j+1,0,t_act)
    cum_econ = economies - dpa
    VANi = va - dpa
    for j in range(1,Nccv):
        econ_combustible = econ_combustible*(1+inf_comb)
        interet = balance*inte
        cap_paye = paiement-interet
        balance  = balance - cap_paye
        costs = paiement + maintenance[j]*(1+inf_ord)**(j)    # flux négatifs non actualisés
        gains = econ_combustible  + deduct_impot*interet                                                 # flux positifs non actualisés
        economies = gains - costs
        va = present_value(economies,j+1,0,t_act)
        cum_econ = cum_econ + economies
        VANi = VANi + va
        if abs(balance) < 1:
            paiement = 0
            balance = 0
    vrn = present_value(vr,Nccv,0,t_act)
    VAN = VANi + vrn
    return VAN

def calcul_miir(rr,fr,dpa,paiementi,econ_combi,vr,mi,maintenance):

    economies = np.zeros(Nccv)
    num = 0

    den = dpa
    paiement = paiementi
    j = 0
    econ_combustible = econ_combi
    interet = mi*inte
    cap_paye = paiement-interet
    balance  = mi - cap_paye
    costs = paiement
    gains = econ_combustible + deduct_impot*interet
    economies[j] = gains - costs
    if economies[j] > 0:
        num = num + future_value(economies[j],N-(j+1),rr)
    else:
        den = den - present_value(economies[j],j+1,0,fr)
    for j in range(1,Nccv):
        econ_combustible = econ_combustible*(1+inf_comb)
        interet = balance*inte
        cap_paye = paiement-interet
        balance  = balance - cap_paye
        costs = paiement + maintenance[j]*(1+inf_ord)**(j)    # flux négatifs non actualisés
        gains = econ_combustible  + deduct_impot*interet                                                 # flux positifs non actualisés
        economies[j] = gains - costs
        if j == Nccv-1:
            economies[j] = economies[j] + vr
        if abs(balance) < 1:
            paiement = 0
            balance = 0
        if economies[j] > 0:
            num = num + future_value(economies[j],Nccv-(j+1),rr)
        else:
            den = den - present_value(economies[j],j+1,0,fr)
    miir = (num/den)**(1/Nccv) - 1
    return miir,economies


# ajout de la consommation électrique de la pompe en suposant qu'elles fonctionnent toujuors
#
# analyse de cycle de vie (relative)
# Cout dela PAC
#
dpaa = Investa*mise_depart
mia =  Investa - dpaa
paiementia = mia/pwf(Npret,0,inte)
econ_caria = cout_base - cout_caria
VANa = calcul_npv(disc,dpaa,paiementia,econ_caria,vra,mia,maintenancea)


#
# 2eme cas
#
dpab = Investb*mise_depart
mib =  Investb - dpab
paiementib = mib/pwf(Npret,0,inte)
#
econ_carib = cout_base - cout_carib
VANb = calcul_npv(disc,dpab,paiementib,econ_carib,vrb,mib,maintenanceb)


def calcul_tri(tt,dpa,paiementi,econ_cari,vrp,mi,glycol_cost=0):
    vann = calcul_npv(tt,dpa,paiementi,econ_cari,vrp,mi,glycol_cost)
    return vann
tria = newton(calcul_tri,.09,args = (dpaa,paiementia,econ_caria,vra,mia,maintenancea))
trib = newton(calcul_tri,.09,args = (dpab,paiementib,econ_carib,vrb,mib,maintenanceb))
VANan = calcul_npv(tria,dpaa,paiementia,econ_caria,vra,mia,maintenancea)
VANbn = calcul_npv(trib,dpab,paiementib,econ_carib,vrb,mib,maintenanceb)
print(tria,trib)
print(VANan,VANbn)
exit()
fr = 0.06
rr = 0.06
miira1,eca = calcul_miir(tria,tria,dpaa,paiementia,econ_caria,vra,mia,maintenancea)
profits = np.maximum(eca,0)
pertes = np.maximum(-eca,0)
miira2 = mirr(dpaa,profits,pertes,tria,tria)
print(tria,miira1,miira2)
miirb1,ecb = calcul_miir(trib,trib,dpab,paiementib,econ_carib,vrb,mib,maintenanceb)
profits = np.maximum(ecb,0)
pertes = np.maximum(-ecb,0)
miirb2 = mirr(dpab,profits,pertes,trib,trib)
print(trib,miirb1,miirb2)

