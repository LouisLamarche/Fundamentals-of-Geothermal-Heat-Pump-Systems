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
To = 9.0
COP_h = 3.4
COP_c = 4.8
ks = 3.1
alj = 0.09
alhr = alj/24.0
qh_ch_bat = 11000.0 # pulse horaire en Watts batiment
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
qm_ch = max(qm)
qh_ch = qh_ch_bat*(COP_h-1)/COP_h

#
# summation of the first  12 premiers pulses
#
# 1er
#
# question b
#     total consumtion  W-jr/m
Cons = sum(qm*hrm)
#
#
qa = Cons/8760.0      # mean
#
# add monthly  pulse
#

#
# summation of the first  12 premiers pulses
#
# 1er
#
# question b
#     total consumtion  W-jr/m
Cons = sum(qm*hrm)
#
#
qa = Cons/8760.0      # mean
#
# add monthly  pulse
#
# temps
#
nannees  = 20
ta = nannees*8760.0
tm = ta + 730.0
tf = tm + 4
Fof = alhr*tf/rb**2
Fo1 = alhr*(tf-ta)/rb**2
Fo2 = alhr*(tf-tm)/rb**2
#
# debut du design
# Calcul des Tf
#
#Tfi_ch = Tfo_ch - qh_ch/(mp*Cp)
Tf_ch = 0
#
# Calul des résistances du sol
#
Gf = G_function(Fof)
G1 = G_function(Fo1)
G2 = G_function(Fo2)
Rsa = (Gf-G1)/ks
Rsm = (G1-G2)/ks
Rsh = (G2)/ks
# calcul de la longeur Chauffage
Tp = 0.0
Li = 300.0
ok = False
compt = 0
delta = 0.1
nx = 2
ny = 1
nt = nx*ny
b = 6
while not ok:
    Lch1 = (qa*Rsa)/(To - Tf_ch - Tp)
    Lch2 = (qm_ch*Rsm)/(To - Tf_ch - Tp)
    Lch3 = (qh_ch*Rsh)/(To - Tf_ch- Tp)
    Lch4 = (qh_ch*Rpb)/(To - Tf_ch - Tp)
    L = Lch1 + Lch2 + Lch3 + Lch4
    # calcul de la longeur Climatisation
    H = L/nt
    print ('length per borehole is ' ,'%.1f' % H ,  ' m')
    qp = qa/L
    Tpn = Tp_Bernier(nx,ny,b,qp,alj,ks,n_annees=20)
    print ('Penalty temperature is ','%.2f' % Tpn, ' C')
    if (abs(L-Li)) < delta:
        ok = True
    else:
        compt = compt + 1
        if compt > 1000:
            print ('divergence')
            ok = True
        Tp = Tpn
        Li = L
Ls = '{:6.2f}'.format(L)   # chaine de caractere avec 2 decimales
print ('la longueur est '+ Ls + ' m')









