#
# exemple solaire 11.4
#
# Exemple de calcul du cout de cycle de vie ( analyse économique)
#
#

from sys import *
from solar_mod import *
from finance_mod import *
from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import newton
from xlwt import *

jm = array([17,47,75,105,135,162,198,228,258,288,318,344])# jour type de chaque mois
nm = array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois

Invest = 38000.0    # cout total
financement = 0.9   # pourcentage de finacement
dpa = Invest*(1.0 - financement)        # paiement initial
mi = Invest - dpa         # montant à payer
vr = 2000.0            # valeur résiduelle
t = 0.08            # taux d'actualisation
inte = 0.04         # intéret sur pret
N = 15              # nombre d,annéees  pour le calcul de la VAN
Npret = 10          # nombre d'annéees du pret
inf_ass = 0.03      # inflation sur assurances
inf_tax = 0.03      # inflation sur taxes
inf_comb = 0.05      # inflatio sur combustible
extra_assi = 0.0  # cout extra assurance
extra_taxi = 0.0  # cout extra tax
cout_cari = 8000.0  # cout de carburant initial
pourc_econ = 0.4   # fraction solaire
deduct_impot = 0.3 # deduction fiscale sur les taxe fonciere ainsi que l'interet
econ_combi = cout_cari*pourc_econ           # economie combustible initial
paiementi = mi/calcul_fva(Npret,0,inte)     # paiment annuel
paiement = paiementi
maintenance = np.zeros(N)
maintenance[4:N:5] = 500.0
#
# initialisation des vecteurs
#
econ_combustible = zeros(N)
va = zeros(N)
extra_taxv = zeros(N)
extra_assv = zeros(N)
balance = zeros(N)
econ_solaire = zeros(N)
cum_econ = zeros(N)     # cumulatif deseconomies  solaires
VANi = zeros(N)
deduct_impv = zeros(N)
paiementv = zeros(N)

#
# analyse de la première année
#
j = 0
cout_combustible = cout_cari
econ_combustible[j] = econ_combi
interet = mi*inte
cap_paye = paiement-interet
balance[j]  = mi - cap_paye
deduct_impv[j] = deduct_impot*(interet+extra_taxi)
totaux_dep_solaire = paiement +  extra_assi + extra_taxi
totaux_gain_solaire = econ_combustible[j] + deduct_impv[j]
econ_solaire[j] = totaux_gain_solaire - totaux_dep_solaire
va[j] = calcul_va(econ_solaire[j],j+1,0,t)
extra_taxv[j] = extra_taxi
extra_assv[j] = extra_assi
cum_gain = totaux_gain_solaire
cum_econ[j] = econ_solaire[j] - dpa
VANi[j] = va[j] - dpa
for j in range(1,N):
    extra_tax = extra_taxv[j-1]*(1+inf_tax)
    extra_ass = extra_assv[j-1]*(1+inf_tax)
    cout_combustible = cout_combustible*(1+inf_comb)
    econ_combustible[j] = cout_combustible*pourc_econ
    interet = balance[j-1]*inte
    cap_paye = paiement-interet
    balance[j]  = balance[j-1] - cap_paye
    paiementv[j] = paiement
    deduct_impv[j] = deduct_impot*(interet + extra_tax)
    mai = maintenance[j]*(1+inf_tax)**(j)
    totaux_dep_solaire = paiement + extra_tax + extra_ass + mai
    totaux_gain_solaire = econ_combustible[j] + deduct_impv[j]
    econ_solaire[j] = totaux_gain_solaire - totaux_dep_solaire
    va[j] = calcul_va(econ_solaire[j],j+1,0,t)
    extra_taxv[j] = extra_tax
    extra_assv[j] = extra_ass
    cum_gain = cum_gain + totaux_gain_solaire
    cum_econ[j] = cum_econ[j-1] + econ_solaire[j]
    VANi[j] = VANi[j-1] + va[j];
    if abs(balance[j]) < 1:
        paiement = 0
        balance[j] = 0
vrn = calcul_va(vr,N,0,t)
VANi[N-1] = VANi[N-1]
VAN = VANi[N-1] + vrn
VAN2 = -paiementi*calcul_fva(Npret,0,t) - extra_taxi*calcul_fva(N,inf_tax,t) - extra_assi*calcul_fva(N,inf_ass,t) \
        + econ_combi*calcul_fva(N,inf_comb,t)  + deduct_impot*extra_taxi*calcul_fva(N,inf_tax,t) + deduct_impot*mi*calcul_fva_int(Npret,N,inte,t) - dpa + vrn \
        - maintenance[4]*(1+inf_tax)**(4)/(1+t)**5 - maintenance[9]*(1+inf_tax)**(9)/(1+t)**10 \
        - maintenance[14]*(1+inf_tax)**(14)/(1+t)**15 \

print ('La valeur actualisée nette est de = ',VAN)
print ('La valeur actualisée nette est de = ',VAN2)
exit()
def calcul_tri(tt):
    vrnn = calcul_va(vr,N,0,tt)
    vann = -paiementi*calcul_fva(Npret,0,tt) - extra_taxi*calcul_fva(N,inf_tax,tt) - extra_assi*calcul_fva(N,inf_ass,tt) \
        + econ_combi*calcul_fva(N,inf_comb,tt) + deduct_impot*extra_taxi*calcul_fva(N,inf_tax,tt)+ deduct_impot*mi*calcul_fva_int(Npret,N,inte,tt) - dpa + vrnn
    return vann
tri = newton(calcul_tri,0.1)
print ('Le tri est de = ',tri)
# verification
t3 = tri
vrn3= calcul_va(vr,N,0,t3)
VAN3 = -paiementi*calcul_fva(Npret,0,t3) - extra_taxi*calcul_fva(N,inf_tax,t3) - extra_assi*calcul_fva(N,inf_ass,t3) \
        + econ_combi*calcul_fva(N,inf_comb,t3)  + deduct_impot*extra_taxi*calcul_fva(N,inf_tax,t3) + deduct_impot*mi*calcul_fva_int(Npret,N,inte,t3) - dpa + vrn3
print ('La VAN pou t = ' + str(t3) + ' est de = ',VAN3)
#
# calcul des temps de retour sur investissment
#
x = arange(0,N)
Inv = Invest*ones(N)
#
#
extra = extra_taxv + extra_assv
flag_excel = 0
if flag_excel:
    classeur = Workbook()
    ic = 0
    feuille1 = classeur.add_sheet('tableau  ')
    for i in range(0,11):         # On agrandit la largeur des colonnes
        col = feuille1.col(i)
        col.width = 256 * 15
    style0 = easyxf('font :bold on')
    style1 = easyxf("pattern: pattern solid, fore_color light_green; font: color black; align: horiz center")
    fmt = '0.0'
    style1.num_format_str = fmt
    # debut de l'écriture
    feuille1.write(ic,0,'Annee  ',style0)
    feuille1.write(ic,1,'Economie combustibles',style0)
    feuille1.write(ic,2,'Paiment annuel',style0)
    feuille1.write(ic,3,'Extra Ass',style0)
    feuille1.write(ic,4,'Extra Taxe',style0)
    feuille1.write(ic,5,'Deduct fiscales',style0)
    feuille1.write(ic,6,'Economies',style0)
    feuille1.write(ic,7,'Economies actualisees',style0)
    feuille1.write(ic,8,'Economies cumulatives',style0)
    feuille1.write(ic,9,'Economies cumulatives actualisees',style0)
    feuille1.write(ic,10,'Balance de capital',style0)
    feuille1.write(ic+1,0,0,style1)
    feuille1.write(ic+1,6,-dpa,style1)
    feuille1.write(ic+1,7,-dpa,style1)
    feuille1.write(ic+1,8,-dpa,style1)
    feuille1.write(ic+1,9,-dpa,style1)
    feuille1.write(ic+1,10,mi,style1)
    for i in range (0,N):      #
        feuille1.write(i+ic+2,0,i+1,style1)
        feuille1.write(i+ic+2,1,econ_combustible[i],style1)
        feuille1.write(i+ic+2,2,paiementv[i],style1)
        feuille1.write(i+ic+2,3,extra_assv[i],style1)
        feuille1.write(i+ic+2,4,extra_taxv[i],style1)
        feuille1.write(i+ic+2,5,deduct_impv[i],style1)
        feuille1.write(i+ic+2,6,econ_solaire[i],style1)
        feuille1.write(i+ic+2,7,va[i],style1)
        feuille1.write(i+ic+2,8,cum_econ[i],style1)
        feuille1.write(i+ic+2,9,VANi[i],style1)
        feuille1.write(i+ic+2,10,balance[i],style1)
    classeur.save('Exemple11_4nouv.xls')
#
# calcul des temps de retour sur investissment
#
