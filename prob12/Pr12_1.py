#coding: utf-8
#
#
from finance_md import *
import numpy as np
from matplotlib.pyplot import *

t = 0.08                # taux d'actualisation
inte = 0.04             # intéret sur pret
N = 20                  # nombre d,annéees  pour le calcul de la VAN
Npret = 10              # nombre d'annéees du pret
inf_comb = 0.04         # inflation sur combustible
Cons = 180*1000/3.6     # en KWh
cout_elec = 0.07        # 7 cents par kWh
cout_cari = cout_elec*Cons      # cout de carburant initial
cout_fixe = 12500
cout_metre = 50
mise = 0.2


def calcul_npv(H,pourc_econ,cout_energie):
    Invest = cout_fixe + H*cout_metre
    dpa = mise*Invest            # paiement initial
    mi = Invest - dpa       # montant résiduel à payer
    econ_combi = cout_energie*pourc_econ           # economie combustible initial
    paiementi = mi/pwf(Npret,0,inte)     # paiment annuel
    VAN = -paiementi*pwf(Npret,0,t)  + econ_combi*pwf(N,inf_comb,t)  - dpa
    return VAN

# cas = 1
H1 = 310.0
pourc_econ1 = 0.6
NPV1 = calcul_npv(H1,pourc_econ1,cout_cari )
print ('NPV  (case 1)  = ',NPV1)
# cas = 2
H2 = 495.0
pourc_econ2 = 0.72
NPV2 = calcul_npv(H2,pourc_econ2,cout_cari )
print ('NPV  (case 2)  = ',NPV2)

Invest1 = cout_fixe + H1*cout_metre
dpa1 = mise*Invest1            # paiement initial1
mi1 = Invest1 - dpa1       # montant résiduel à payer
paiementi1 = mi1/pwf(Npret,0,inte)     # paiment annuel
Invest2 = cout_fixe + H2*cout_metre
dpa2 = mise*Invest2            # paiement initial1
mi2 = Invest2 - dpa2       # montant résiduel à payer
paiementi2 = mi2/pwf(Npret,0,inte)     # paiment annuel

x1 = -paiementi1*pwf(Npret,0,t)   - dpa1
x2 = -paiementi2*pwf(Npret,0,t)   - dpa2
cout_elec = (x1 - x2)/((pourc_econ2 - pourc_econ1)*pwf(N,inf_comb,t)*Cons)
print('Electricity cost = ',cout_elec)
#
# verification
cout_carib =cout_elec*Cons
NPV1 = calcul_npv(H1,pourc_econ1,cout_carib)
print ('NPV  (case 1)  = ',NPV1)
NPV2 = calcul_npv(H2,pourc_econ2,cout_carib)
print ('NPV  (case 2)  = ',NPV2)




