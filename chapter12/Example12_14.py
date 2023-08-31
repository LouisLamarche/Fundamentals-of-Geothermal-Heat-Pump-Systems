#
#
from finance_md import *
import numpy as np
from matplotlib.pyplot import *
from scipy.optimize import newton
from xlwt import *

Invest = 30000.0        # cout total
dpa = 0.1*Invest        # paiement initial
mi = Invest - dpa       # montant résiduel à payer
vr = 1000                # valeur résiduelle
t = 0.08                # taux d'actualisation
inte = 0.03             # intéret sur pret
N = 20                  # nombre d,annéees  pour le calcul de la VAN
Npret = 10              # nombre d'annéees du pret
inf_tax = 0.03          # inflation sur taxes
inf_comb = 0.045        # inflation sur combustible
extra_taxi = 100.0      # cout extra tax
cout_cari = 3500.0      # cout de carburant initial
deduct_imp = 0.00       # deduction d'impot sur l'achat et l'augmentation de l'assurance
pourc_econ = 0.60
econ_combi = cout_cari*pourc_econ           # economie combustible initial
paiementi = mi/pwf(Npret,0,inte)     # paiment annuel
paiement = paiementi
#
econ_combustible = np.zeros(N)
cap_paye = np.zeros(N)
paiementv = np.zeros(N)
extra_taxv = np.zeros(N)
balance = np.zeros(N)         # solde du pret
gains = np.zeros(N)           # gains moins depenses
costs = np.zeros(N)           # gains moins depenses
economies = np.zeros(N)           # gains moins depenses
cum_econ = np.zeros(N)     # cumulatif des economies
va = np.zeros(N)           # valeur actualisee des economies
VANi = np.zeros(N)         # valeur actualisee nette cumulative
#
fr = 0.05
rr = 0.07
maintenance = np.zeros(N)
maintenance[9] = 300
#
# initialisation des variables
num = 0
den = dpa
j = 0
econ_combustible[j] = econ_combi
interet = mi*inte
cap_paye[j] = paiement-interet
balance[j]  = mi - cap_paye[j]
extra_taxv[j] = extra_taxi
costs[j] = paiement + extra_taxv[j]
gains[j] = econ_combustible[j]
economies[j] =gains[j] - costs[j]
va[j] = present_value(economies[j],j+1,0,t)
cum_econ[j] = economies[j] - dpa
VANi[j] = va[j] - dpa
paiementv[j] = paiement
if economies[j] > 0:
    num = num + future_value(economies[j],N-(j+1),rr)
else:
    den = den - present_value(economies[j],j+1,0,fr)
for j in range(1,N):
    extra_taxv[j] = extra_taxv[j-1]*(1+inf_tax)
    econ_combustible[j] = econ_combustible[j-1]*(1+inf_comb)
    interet = balance[j-1]*inte
    cap_paye[j] = paiement-interet
    balance[j]  = balance[j-1] - cap_paye[j]
    costs[j] = paiement + extra_taxv[j] + maintenance[j]*(1+inf_tax)**(j)    # flux négatifs non actualisés
    gains[j] = econ_combustible[j]                                                   # flux positifs non actualisés
    economies[j] = gains[j] - costs[j]
    va[j] = present_value(economies[j],j+1,0,t)
    cum_econ[j] = cum_econ[j-1] + economies[j]
    VANi[j] = VANi[j-1] + va[j]
    paiementv[j] = paiement
    if abs(balance[j]) < 1:
        paiement = 0
        balance[j] = 0
    if economies[j] > 0:
        num = num + future_value(economies[j],N-(j+1),rr)
    else:
        den = den - present_value(economies[j],j+1,0,fr)

vrn = present_value(vr,N,0,t)
VANi[N-1] = VANi[N-1] + vrn
npv1 = VANi[N-1]
npv2 = -paiementi*pwf(Npret,0,t) - extra_taxi*pwf(N,inf_tax,t) \
        + econ_combi*pwf(N,inf_comb,t) - dpa + vrn - maintenance[9]*(1+inf_tax)**(9)/(1+t)**10
print (' The  NPV is = ',npv1)
print (' The  NPV is = ',npv2)

def calcul_tri(tt):
    vrnn = present_value(vr,N,0,tt)
    vann = -paiementi*pwf(Npret,0,tt) - extra_taxi*pwf(N,inf_tax,tt) \
        + econ_combi*pwf(N,inf_comb,tt) - dpa + vrnn - maintenance[9]*(1+inf_tax)**(9)/(1+tt)**10
    return vann
tri = newton(calcul_tri,t)
print ('the iir is = ',tri)
# verification
t = tri
vrn = present_value(vr,N,0,t)
npv3 = -paiementi*pwf(Npret,0,t) - extra_taxi*pwf(N,inf_tax,t) \
        + econ_combi*pwf(N,inf_comb,t) - dpa + vrn - maintenance[9]*(1+inf_tax)**(9)/(1+t)**10
print('THe NPV for t = ' + str(tri) +' is = ' + str(npv3))
flag_excel = 0
if flag_excel:
    classeur = Workbook()
    ic = 0
    feuille1 = classeur.add_sheet('tableau  ')
    for i in range(0,9):         # On agrandit la largeur des colonnes
        col = feuille1.col(i)
        col.width = 256 * 15
    fmt = '0'
    style0 = easyxf('font :bold on')
    style1 = easyxf("pattern: pattern solid, fore_color light_green; font: color black; align: horiz center")
    style1.num_format_str = fmt
    # debut de l'écriture
    feuille1.write(ic,0,'Year  ',style0)
    feuille1.write(ic,1,'Energy gain',style0)
    feuille1.write(ic,2,'Annual payment',style0)
    feuille1.write(ic,3,'Cash flow',style0)
    feuille1.write(ic,4,'discounted profit',style0)
    feuille1.write(ic,5,'Cumulative cash flow ',style0)
    feuille1.write(ic,6,'Discouted cumulative cash flow ',style0)
    feuille1.write(ic,7,'Balance',style0)
    feuille1.write(ic+1,0,0,style1)
    feuille1.write(ic+1,6,-dpa,style1)
    feuille1.write(ic+1,7,-dpa,style1)
    feuille1.write(ic+1,8,mi,style1)
    for i in range (0,N):      #
        feuille1.write(i+ic+2,0,i+1,style1)
        feuille1.write(i+ic+2,1,econ_combustible[i],style1)
        feuille1.write(i+ic+2,2,paiementv[i],style1)
        feuille1.write(i+ic+2,3,extra_taxv[i],style1)
        feuille1.write(i+ic+2,4,economies[i],style1)
        feuille1.write(i+ic+2,5,va[i],style1)
        feuille1.write(i+ic+2,6,cum_econ[i],style1)
        feuille1.write(i+ic+2,7,VANi[i],style1)
        feuille1.write(i+ic+2,8,balance[i],style1)
    classeur.save('Exemple12_14.xls')
profits = np.maximum(economies,0)
pertes = np.maximum(-economies,0)
miira = (num/den)**(1/N) - 1
miirb = mirr(dpa,economies,rr,fr)
print ('miir is = ',miira,miirb)
x = np.arange(1,N+1)

p1 = plot(x,cum_econ,'k-x',x,VANi,'k-o',x,balance,'k-+',markersize=11)
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
rr = 0.05
grid(True,which='both')
fts = 16
ftst = 14
legend(('Cumulative cash flow','Cumulative discounted cash flow','Mortgage balance'),fontsize = ftst)
xlabel('Year',fontsize = fts)
#ylabel(r'LCC',fontsize = fts)
xticks(fontsize=ftst)
yticks(fontsize=ftst)
show()


