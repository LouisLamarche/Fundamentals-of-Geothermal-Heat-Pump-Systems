#coding: utf-8
#
# Exemple 10.2 (Conception de puits horizontaux)
#
from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
# vecteurs d'heures
hrm = array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = array([ 744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])

def Calcul_R_horizontal_gauche(X1,zt1,zt2,ks):
    X2 = (zt2-zt1)*X1
    X3 = 2*zt1*X1
    X4 = (zt1+zt2)*X1
    R1 = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    X3 = 2*zt2*X1
    R2 = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    R = (R1+R2)/2.0
    return R

def Calcul_R_horizontal_droite(X1,zt1,xt2,ks):
    X2 = xt2*X1
    X3 = 2*zt1*X1
    dd = sqrt((2*zt1)**2 + xt2**2)
    X4 = dd*X1
    R = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    return R

# sol
ks = 1.2
alj =0.05
alhr = alj/24
COP_ch = 3.0
COP_cl = 4.0

#
# conduite
#
Rpp = 0.08
d2 = 0.05        # diametre exterieur
rp = d2/2.0

#
#
# fluide
Patm = 101.325*1000.0
cas = 2
if cas == 1:
    fluide = 'Water'
elif cas ==2:
    fluide  = 'INCOMP::APG-30%'   # propylene - glycol 30 %
Trefk = 273 + 30
rho = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Pr = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)

mpr = 0.054/1000.0*rho  # 0.054 l/s par kW
Tfo_ch = -2.0    # valeur de design
Tfo_cl = 30.0   # valeur de design
q_bat_ch = 8000.0
q_bat_cl = -10000.0
CAPCH = q_bat_ch/1000
CAPCL = -q_bat_cl/1000
mp_ch = mpr*CAPCH        # debit
mp_cl = mpr*CAPCL        # debit
#
# charges
#
PLFmcl = 0.25
PLFmch = 0.1
EFLH_ch = 400
EFLH_cl = 800
qh_ch = q_bat_ch*(COP_ch-1)/COP_ch
qh_cl = q_bat_cl*(COP_cl+1)/COP_cl
qm_ch = qh_ch*PLFmch
qm_cl = qh_cl*PLFmcl
qa = (qh_ch*EFLH_ch + qh_cl*EFLH_cl)/8760.0
#
# temps
#
ta = 10.0*8760.0
tm = ta + 730.0
tf = tm + 4.0
#
# debut du design0
# Calcul des Tf
#
Tfi_ch = Tfo_ch - qh_ch/(mp_ch*Cp)
Tfi_cl = Tfo_cl - qh_cl/(mp_cl*Cp)
Tf_ch = (Tfi_ch+Tfo_ch)/2.0
Tf_cl = (Tfi_cl+Tfo_cl)/2.0
#
# Calul des résistances du sol
#
z = 1.5             # profondeur
zt = z/rp
Xf = rp/(2*sqrt(alhr*tf))
X1 = rp/(2*sqrt(alhr*(tf-ta)))
X2 = rp/(2*sqrt(alhr*(tf-tm)))
#
# cas figure gauche
#
z1 = z - 0.5
zm = (z + z1)/2.0
zt1 = z1/rp
Rhf = Calcul_R_horizontal_gauche(Xf,zt1,zt,ks)
Rh1 = Calcul_R_horizontal_gauche(X1,zt1,zt,ks)
Rh2 = Calcul_R_horizontal_gauche(X2,zt1,zt,ks)
Rsa = Rhf - Rh1
Rsm = Rh1 - Rh2
Rsh = Rh2
print('Rsa (gauche) =',Rsa)
print('Rsm (gauche) =',Rsm)
print('Rsh (gauche) =',Rsh)

#
# calcul de temperature de la Terre
To = 11.0
tshift = 5     # Valeur pour Montreal ( en jour)
amp = 14      # amplitude de variation de la temperature dans l'annéee
tch = 5         #  jour le plus froid
tcl = 180         #  jour le plus chaud
dTchg,dTmaxg = Delta_T_earth(zm,tch,tshift,alj,amp)  # diffusivite en m2/jr
dTclg,dTmaxg = Delta_T_earth(zm,tcl,tshift,alj,amp)
Toch = To + dTchg
Tocl = To + dTclg
#
# calcul de la longeur Chauffage
# 1ere iteration
Lch1 = (qa*Rsa)/(Toch - Tf_ch )
Lch2 = (qm_ch*Rsm)/(Toch - Tf_ch )
Lch3 = (qh_ch*Rsh)/(Toch - Tf_ch )
Lch4 = (qh_ch*Rpp)/(Toch - Tf_ch )
L_ch = Lch1 + Lch2 + Lch3 + Lch4
# calcul de la longeur Climatisation
Lcl1 = (qa*Rsa)/(Tocl - Tf_cl )
Lcl2 = (qm_cl*Rsm)/(Tocl - Tf_cl )
Lcl3 = (qh_cl*Rsh)/(Tocl - Tf_cl)
Lcl4 = (qh_cl*Rpp)/(Tocl - Tf_cl)
L_cl = Lcl1 + Lcl2 + Lcl3 + Lcl4
L = max(L_ch,L_cl)
print('L chauffage =',L_ch/2,'L climatisation = ',L_cl/2)
print ('la longueur circuit gauche est '+ str(L/2) + ' m')

#
# figure droite
#
zm = z
xt2 = 0.5/rp
Rhf = Calcul_R_horizontal_droite(Xf,zt,xt2,ks)
Rh1 = Calcul_R_horizontal_droite(X1,zt,xt2,ks)
Rh2 = Calcul_R_horizontal_droite(X2,zt,xt2,ks)
Rsa = Rhf - Rh1
Rsm = Rh1 - Rh2
Rsh = Rh2
print('Rsa (droite) =',Rsa)
print('Rsm (droite) =',Rsm)
print('Rsh (droite) =',Rsh)
#
# calcul de temperature de la Terre
To = 11.0
tshift = 5     # Valeur pour Montreal ( en jour)
amp = 14      # amplitude de variation de la temperature dans l'annéee
tch = 5         #  jour le plus froid
tcl = 180         #  jour le plus chaud
dTchd,dTmaxd = Delta_T_earth(zm,tch,tshift,alj,amp)  # diffusivite en m2/jr
dTcld,dTmaxd = Delta_T_earth(zm,tcl,tshift,alj,amp)
Toch = To + dTchd
Tocl = To + dTcld
#
# calcul de la longeur Chauffage
# 1ere iteration
Lch1 = (qa*Rsa)/(Toch - Tf_ch )
Lch2 = (qm_ch*Rsm)/(Toch - Tf_ch )
Lch3 = (qh_ch*Rsh)/(Toch - Tf_ch )
Lch4 = (qh_ch*Rpp)/(Toch - Tf_ch )
L_ch = Lch1 + Lch2 + Lch3 + Lch4
# calcul de la longeur Climatisation
Lcl1 = (qa*Rsa)/(Tocl - Tf_cl )
Lcl2 = (qm_cl*Rsm)/(Tocl - Tf_cl )
Lcl3 = (qh_cl*Rsh)/(Tocl - Tf_cl)
Lcl4 = (qh_cl*Rpp)/(Tocl - Tf_cl)
L_cl = Lcl1 + Lcl2 + Lcl3 + Lcl4
L = max(L_ch,L_cl)
print('L chauffage =',L_ch/2,'L climatisation = ',L_cl/2)
print ('la longueur circuit droite est '+ str(L/2) + ' m')

