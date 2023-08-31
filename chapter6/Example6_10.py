#
# Exemple 6.10 dimensionnemnt approche suédoise
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
#
# data
#
cas = 'a'
# soil
als = 0.09 # m2/jr
alhr = als/24.0
ks = 2.5
To = 10.0
COP_ch = 3.0
COP_cl = 5.0
CAP_ch = 12000.0*3
CAP_cl = 12000.0*3
# well
rb = 0.15/2.0
kg = 1.7
do = 0.033
di = 0.027
ro = do/2.0
ri = di/2.0
xc = rb/3
kp = 0.4
zo = 4.0
# fluid
rhof = 1030
mp1 = 0.6/1000.0*rhof       # flow rate per borehole
Cpf = 3800.0
muf = 0.004
kf = 0.44
alf = kf/(rhof*Cpf)
nuf = muf/rhof
Pr = nuf/alf
#
# design data
Tfo_ch = 0.0    # valeur de design
Tfo_cl = 25.0   # valeur de design
n_years = 10
n_rings = 3
nh = 4          # hourly block
#
# field configuration
nx = 3
ny = 2
b = 6.0
b = 6.0
nt = nx*ny
mp = mp1*nt
# resistance
Re = 4*mp1/(pi*di*muf)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hf = (Nud*kf)/(di)
rconv = 1/(pi*di*hf)
rcond = np.log(do/di)/(2*pi*kp)
Rp = rcond + rconv
Rb1 = Rb_Paul(kg,rb,ro,'b') + Rp/2.0
J = 10
z = np.array([xc,-xc])  # pipes coordinates
Rb2,Ra2 = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
#mois_ini = 5 # june
if cas == 'a':
    Rb = Rb1
    mois_ini = 0 # january
elif cas == 'b':
    Rb = Rb1
    mois_ini = 5 # june
elif cas == 'c':
    Rb = Rb2
    mois_ini = 0 # january
else:
    Rb = Rb2
    mois_ini = 5 # june

hrm = 730
A = np.array([[2.17,0.00],[2.07, 0.00],[1.750,0.000],[1.39,0.00],[0.9,-1.2],[0.00,-1.6],[0.00,-2.00],[0.00,-1.80],[0.85,-0.80],[1.22,-0.40],[1.64,0.00],[2.02,0.00]])
A = 3*A     # all loads are multiplied by 3
E_chau = np.sum(A[:,0])  # MWh
q_chau = A[:,0]*1.0e6/730.0  # W
E_clim = np.sum(A[:,1])  # MWh
q_clim = A[:,1]*1.0e6/730.0    # W
q_sol = q_chau+q_clim
qm_ch = max(q_sol)
qm_cl = min(q_sol)
qa = (E_chau + E_clim)*1.0e6/8760   # W
qh_ch = CAP_ch*(COP_ch-1)/COP_ch
qh_cl = -CAP_cl*(COP_cl+1)/COP_cl
#

#
# debut du design
# Calcul des Tf
#
Tfi_ch = Tfo_ch - qh_ch/(mp*Cpf)
Tfi_cl = Tfo_cl - qh_cl/(mp*Cpf)
Tf_ch = (Tfi_ch+Tfo_ch)/2.0
Tf_cl = (Tfi_cl+Tfo_cl)/2.0
#
#############################################################################################
#
# methode suedoise
#
# choix initial de la hauteur de puits
nbloc = 4
Hi = 110.0


ok = False
compt = 1
compt_max = 20
delta = 0.5
mois_ch = np.argmax(q_sol)
mois_cl = np.argmin(q_sol)
if qa > 0:
    nyears_ch = n_years
    nyears_cl = 0
    if mois_cl - mois_ini < 0:
        nyears_cl = 1
else:
    nyears_cl = n_years
    nyears_ch = 0
    if mois_ch - mois_ini < 0:
        nyears_ch = 1
ntot_ch = int(nyears_ch*12 + mois_ch + 1 - mois_ini)  + 1
ntot_cl = int(nyears_cl*12 + mois_cl + 1 - mois_ini)  + 1
tf_ch = (ntot_ch-1)*730.0 + nbloc  # temps final en heures
tf_cl = (ntot_cl-1)*730.0 + nbloc  # temps final en heures
tf = max(tf_ch,tf_cl)

while not ok :
    rr = rb/Hi
    zob = zo/Hi
    ts = Hi**2/(9*alhr)         # temps caractéristique Eskilson en heures
    Bt = b/Hi
    #
    # génération du champ de capteur rectangulaire 3x2 , configuration initiale
    #
    zt = np.zeros([nt,2])
    k = 0
    for i1  in range(0,nx):
        x = i1*Bt
        for i2 in range(0,ny):
            y = i2*Bt
            zt[k] = [x,y]
            k = k+1
    #
    #############################################################################################
    #############################################################################################
    #
    # calcul de la fonction g analytique associée au champ de capteur choisi
    #
    #
    ng = 80
    t1 = 1.0/ts         # temps minimum 1 heure
    t2 = (1.5*tf)/ts    # temps maximum 1.5 fois tfinal
    ttt = np.logspace(np.log10(t1),np.log10(t2),ng)      # On calcule ng valeurs de la fonction g qu'on interpole par la suite
    ntt = ng+1
    tt = np.zeros(ntt)
    g = np.zeros(ntt)
    tt[1:ntt] = ttt[0:ng]
    for i in range (0,ntt):
        g[i] = compute_g_function(zt,tt[i],rbb = rr,zob = zob)    # SI on n'ecrit pas D/H , la fonction prend 0.04 par défaut ( Comme Eskilson)
    #
    #############################################################################################
    #############################################################################################
    #
    # sommation des pulses mensuels
    #
    q1 = q_sol[mois_ini]
    Foo = tf_ch/ts
    SCM  = q1*np.interp(Foo,tt,g)
    ti = 0
    for ia in range(1,ntot_ch-1):
       nm = (ia + mois_ini) % 12   # On trouve le mois dans l'annee
       ti = ti + 730.0
       Foo = (tf_ch-ti)/ts
       SCM = SCM + (q_sol[nm]-q1)*np.interp(Foo,tt,g)
       q1 = q_sol[nm]
    Foo = nbloc/ts
    #
    # ajout du pulse horaire
    #
    SCT = SCM + (qh_ch-q1)*np.interp(Foo,tt,g)
    #
    # ajout de l'effet de la résistance de puits
    #
    num = SCT/(2*pi*ks)+qh_ch*Rb
    Lch_Es = num/(To-Tf_ch)
    # clim
    q1 = q_sol[mois_ini]
    Foo = tf_cl/ts
    SCT  = q1*np.interp(Foo,tt,g)
    ti = 0
    for ia in range(1,ntot_cl-1):
       nm = (ia + mois_ini) % 12   # On trouve l mois dans l'annee
       ti = ti + 730.0
       Foo = (tf_cl-ti)/ts
       SCT = SCT+ (q_sol[nm]-q1)*np.interp(Foo,tt,g)
       q1 = q_sol[nm]
    Foo = nbloc/ts
    # ajout du pulse horaire
    SCT = SCT + (qh_cl-q1)*np.interp(Foo,tt,g)
    # ajout de l'effet de la résistance de puits
    num = SCT/(2*pi*ks)+qh_cl*Rb
    Lcl_Es = num/(To-Tf_cl)
    L_Es = max(Lch_Es,Lcl_Es)
    print ('L_heating = ', Lch_Es)
    print ('L_cooling = ', Lcl_Es)
    H = L_Es/nt
    print ('L = ', L_Es)
    print ('H = ', H)
    if abs(H-Hi) < delta:
        ok = True
    else:
        Hi = H
        compt = compt + 1
        if compt> compt_max:
            err =1
            ok = True

