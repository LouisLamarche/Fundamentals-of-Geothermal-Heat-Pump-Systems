#coding: latin-1
#
# Exemple 10.1 (Conception de puits horizontaux)
#
from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
from time import *
# vecteurs d'heures
hrm = array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = array([ 744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])


def Calcul_R_horizontal_4pipes(X1,zt1,zt2,dt2,ks):

    d1 = zt2 - zt1
    X2 = d1*X1
    X3 = dt2*X1
    X4 = sqrt(dt2**2+d1**2)*X1
    X5 = 2*zt2*X1
    X6 = (zt2+zt1)*X1
    X7 = sqrt(dt2**2+(2*zt2)**2)*X1
    X8 = sqrt(dt2**2+(zt1+zt2)**2)*X1
    Rs1 = ((I_function(X1)+I_function(X2)+I_function(X3)+I_function(X4))-(I_function(X5)+I_function(X6)+I_function(X7)+I_function(X8)))/(2*pi*ks)     # mK/W
    X5 = 2*zt1*X1
    X7 = sqrt(dt2**2+(2*zt1)**2)*X1
    Rs2 = ((I_function(X1)+I_function(X2)+I_function(X3)+I_function(X4))-(I_function(X5)+I_function(X6)+I_function(X7)+I_function(X8)))/(2*pi*ks)     # mK/W
    R = (Rs1+Rs2)/2.0
    return R

def Calcul_R_horizontal_2pipes(X1,zt,xt,ks):

    X2 = xt*X1
    X3 = 2*zt*X1
    X4 = sqrt(xt**2+(2*zt)**2)*X1
    R = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    return R



# sol
ks = 1.5
CC = 2.0e6
al = ks/CC
alhr = al*3600
als = alhr*24
To = 10.0
COP_ch = 4.0
COP_cl = 5.0
CAP_ch = 40000.0
CAP_cl = 20000.0
#
# conduite
#
d2 = 0.035          # diametre exterieur
d1 = 0.03         # diametre exterieur
rp = d2/2.0
z1 = 1.75             # profondeur
d = 0.3
kp = 0.4
h = 5000
Rpp = log(d2/d1)/(2*pi*kp) + 1/(pi*d1*h)
#
#
# fluide
Vp = 2
rhof = 1030
Cpf = 3800
mp = Vp*rhof/1000
CCf = mp*Cpf
Tfo_ch = -2.0    # valeur de design
Tfo_cl = 30.0   # valeur de design
#
# charges
#
A = array([[15,0.00],[10, 0.00],[10,0.000],[6,0.00],[2,-5],[0.00,-8],\
        [0.00,-10.00],[0.00,-10.00],[1.5,-3],[5,0.0],[10,0.00],[10,0.00]])
E_chau = sum(A[:,0])*(COP_ch-1)/COP_ch  # MWh
q_chau = A[:,0]*1.0e6/hrm*(COP_ch-1)/COP_ch    # W
E_clim = sum(A[:,1])*(COP_cl+1)/COP_cl  # MWh
q_clim = A[:,1]*1.0e6/hrm*(COP_cl+1)/COP_cl    # W
q_sol = q_chau+q_clim
qm_ch = max(q_sol)
qm_cl = min(q_sol)
qa = (E_chau + E_clim)*1.0e6/8760   # W
qh_ch = CAP_ch*(COP_ch-1)/COP_ch
qh_cl = -CAP_cl*(COP_cl+1)/COP_cl
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
Tfi_ch = Tfo_ch - qh_ch/CCf
Tfi_cl = Tfo_cl - qh_cl/CCf
Tf_ch = (Tfi_ch+Tfo_ch)/2.0
Tf_cl = (Tfi_cl+Tfo_cl)/2.0
#
# Calul des résistances du sol
#

zt1 = z1/rp
xt = d/rp
Xf = rp/(2*sqrt(alhr*tf))
X1 = rp/(2*sqrt(alhr*(tf-ta)))
X2 = rp/(2*sqrt(alhr*(tf-tm)))
# autre facon
Rhf = Calcul_R_horizontal_2pipes(Xf,zt1,xt,ks)
Rh1 = Calcul_R_horizontal_2pipes(X1,zt1,xt,ks)
Rh2 = Calcul_R_horizontal_2pipes(X2,zt1,xt,ks)
Rsa = Rhf - Rh1
Rsm = Rh1 - Rh2
Rsh = Rh2

#
# calcul de temperature de la Terre
tshift = 15     # Valeur pour Montreal ( en jour)
amp = 8
tch = 15         #  jour le plus froid
tcl = int(tshift+365/2)        #  jour le plus chaud
zm = z1
dTch,dTmax = Delta_T_earth(zm,tch,tshift,als,amp)  # diffusivite en m2/jr
dTcl,dTmax = Delta_T_earth(zm,tcl,tshift,als,amp)
flag_tmin = 1
if flag_tmin ==1:
    Toch = To + dTch
    Tocl = To + dTcl
else:
    Toch = To - dTmax
    Tocl = To + dTmax
#
# calcul de la longeur Chauffage
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
L_cl = 0+ Lcl2 + Lcl3 + Lcl4
L = max(L_ch,L_cl)
print ('la longueur de tranchee 2 tuyaux est '+ str(L/2) + ' m')
print ('la longueur de de conduite  2 tuyaux est '+ str(L) + ' m')

#8
# slinky
#
# processus itératif
# hypothese L = L/10 1 tuyau
p = 0.5  # pitch
Dring = 1.0
Rring = Dring/2
Li = L/2
ok = False
compt_max = 10
compt = 1
delta = 1
Fof = alhr*tf/Rring**2
Fo1 = alhr*(tf-ta)/Rring**2
Fo2 = alhr*(tf-tm)/Rring**2
zt = z1/Rring
r1 = 6.5/Rring
r2 = 20/Rring
while not ok:
    Nt = int((Li-Dring)/p)
    pt = p/Rring
    Rb =  zeros([Nt,2])
    rpt = rp/Rring
    for i in range(0,Nt):
        Rb[i] = [i*pt,0]
    long = (Nt-1)*pt
    Gringf = G_function_slinky(Rb,Fof,rpt,zt,pt,1,r1,r2)
    Gring1 = G_function_slinky(Rb,Fo1,rpt,zt,pt,1,r1,r2)
    Gring2 = G_function_slinky(Rb,Fo2,rpt,zt,pt,1,r1,r2)
    Rsa_sl = (Gringf - Gring1)/ks
    Rsm_sl = (Gring1 - Gring2)/ks
    Rsh_sl = Gring2/ks
    Rp_sl = Rpp*pt/(2*pi+2*pt)
    # 1ere iteration
    Lch1 = (qa*Rsa_sl)/(Toch - Tf_ch )
    Lch2 = (qm_ch*Rsm_sl)/(Toch - Tf_ch )
    Lch3 = (qh_ch*Rsh_sl)/(Toch - Tf_ch )
    Lch4 = (qh_ch*Rp_sl)/(Toch - Tf_ch )
    L_chsl = Lch1 + Lch2 + Lch3 + Lch4
    # calcul de la longeur Climatisation
    Lcl1 = (qa*Rsa_sl)/(Tocl - Tf_cl )
    Lcl2 = (qm_cl*Rsm_sl)/(Tocl - Tf_cl )
    Lcl3 = (qh_cl*Rsh_sl)/(Tocl - Tf_cl)
    Lcl4 = (qh_cl*Rp_sl)/(Tocl - Tf_cl)
    L_clsl = Lcl1 + Lcl2 + Lcl3 + Lcl4
    Ln = max(L_chsl,L_clsl)
    if (abs(Ln-Li)) < delta:
        ok = True
    else:
        Li = Ln
        compt = compt +1
        if compt > compt_max:
            print ('erreur')
            ok = True
# longeur tranchée
Lt = Nt*p + Dring
# longeur tuyau
Lp = Nt*(pi*Dring+2*p)+Dring*(pi/2+1);print('Lpipe = ',Lp)
print ('la longueur de tranchee slinky est '+ str(Lt) + ' m')
print ('la longueur de de conduite  slinky est '+ str(Lp) + ' m')
Lglhe = 2640
Hglhe = 319
Ntglhe = (Hglhe - Dring)/p
# longeur tuyau
Lplge = Ntglhe*(pi*Dring+2*p)+Dring*(pi/2+1);print('Lpipe = ',Lp)
print ('la longueur de tranchee slinky est '+ str(Hglhe) + ' m')
print ('la longueur de de conduite  slinky est '+ str(Lplge) + ' m')



