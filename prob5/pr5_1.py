#coding: utf-8
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = np.array([ 744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])
#
#
# data
#
Rpb = 0.12
rb = 0.08
H = 100.0
To = 9.0
COP_h = 3.4
COP_c = 4.8
ks = 3.1
alj = 0.09
alhr = alj/24.0
qh_ch = 11000.0 # pulse horaire en Watts batiment
bloc = 4.0          # bloc horaire
tf = 8760 + hrm[0] + bloc   # temps en heures
A = np.array([[9.,	0],
[7.,	0],
[5.,	0],
[4.,	0],
[2.,	0],
[0.,	0],
[0.,	2],
[0.,	3.5],
[0.,	0.5],
[2.,	0],
[5.,	0],
[6.,	0],
])
#
#
#
qch = A[:,0]*1000.0
qcl = A[:,1]*1000.0
#
# ground loads
#
qsol_ch = qch*(COP_h-1)/COP_h
qsol_cl = -qcl*(COP_c+1)/COP_c
#
qm = (qsol_ch + qsol_cl)
qp = qm/H
qph = qh_ch*(COP_h-1)/COP_h/H
#
# summation of the first  12 premiers pulses
#
# 1er
#
Fof = alhr*tf/rb**2
sum1 = qp[0]*G_function(Fof)/ks
ti = 0
#
# sommation des 11 suivants
#
for i in range(1,12):
    ti = ti+hrm[i-1]
    Foi = alhr*(tf - ti)/rb**2
    sum1 = sum1 +  (qp[i]-qp[i-1])*G_function(Foi)/ks
#
# add january (coldest month)
#
Foi = alhr*(tf - 8760.0)/rb**2
sum2 = sum1+ (qp[0]-qp[11])*G_function(Foi)/ks
#
# add hourly
#
Foi = alhr*bloc/rb**2
DeltaTa = sum2+(qph-qp[0])*G_function(Foi)/ks
Tb = To-DeltaTa
Tf = Tb - qph*Rpb
print('a)')
print ('\t Fluid temperature   a)  ',Tf,' C')
#
# question b
#     total consumtion  W-jr/m
Cons = sum(qp*hrm)
#
#
#
qa = Cons/8760.0      # mean
print('b)')
print ('\t annual pulse is   ',qa, ' W/m')
Fof = alhr*tf/rb**2
sum1b = qa*G_function(Fof)/ks
#
# add monthly  pulse
#
Foi = alhr*(tf - 8760.0)/rb**2
sum2b = sum1b + (qp[0]-qa)*G_function(Foi)/ks
#
# add hourly pulse
#
Foi = alhr*bloc/rb**2
DeltaTb = sum2b +(qph-qp[0])*G_function(Foi)/ks
Tbb = To- DeltaTb
Tfb = Tbb - qph*Rpb               # valeur de Rb vient de la question 1
print ('\t Fluid temperature   b)  ',Tfb,' C')
exit()
#

#
# valeur souhaitée
#
Tfn = -5.0
#
# Calcul de SCT
#
SCTa = DeltaTa*H   #  Eq 4.26  en utilisant 14 pulses
SCTb = DeltaTb*H   #  Eq 4.26  en utilisant 3 pulses
#
# Formule 4.26 des notes
Ha = (SCTa + qph*H*Rb)/(To-Tfn)
print('d)')
print ('\t Hauteur nécessaire en utilisant 14 pulses (Tb trouvé en a) ' , Ha, ' m')
Hb = (SCTb + qph*H*Rb)/(To-Tfn)
print (' \t Hauteur nécessaire en utilisant 3 pulses (Tb trouvé en b) ' , Hb, ' m')
#
# Vérification pour le b)
#
print('')
print('Vérification')
print('')
qp = qm/Hb
qph = qh_ch*(COP_ch-1)/COP_ch/Hb
#   calcul de la consommation  totale
co = qp*hrm      # consommation mensuelle en W-hr/m pour chaque mois
Cons = sum(co) # consommation totale sur l'annéee en W-jr/m
qa = Cons/8760.0      # puissance moyenne en W/m
Fof = alhr*tf/rb**2
sum1b = qa*g_function(Fof)/ks
#
# ajout du pulse mensuel
#
Foi = alhr*(tf - 8760.0)/rb**2
sum2b = sum1b + (qp[0]-qa)*g_function(Foi)/ks
#
# ajout du pulse horaire
#
Foi = alhr*bloc/rb**2
DeltaTc = sum2b +(qph-qp[0])*g_function(Foi)/ks
Tbc = To- DeltaTc
print (('\t La  temperature de puits finale est   en c)  ' + str(Tbc) + ' C'))
Tf = Tbc - qph*Rb               # valeur de Rb vient de la question 1
print (('\t La  temperature de fluide finale en c) est    ' + str(Tf) + ' C'))
#