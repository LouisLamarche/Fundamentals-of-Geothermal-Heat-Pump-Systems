# Exemple 6.9 written by Louis lamarche 3 october 2017
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
#
# data
#
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
Rb = Rb_Paul(kg,rb,ro,'b') + Rp/2.0
#
# loads MWh
#
hrm = 730
A = np.array([[2.17,0.00],[2.07, 0.00],[1.750,0.000],[1.39,0.00],[0.9,-1.2],[0.00,-1.6],[0.00,-2.00],[0.00,-1.80],[0.85,-0.80],[1.22,-0.40],[1.64,0.00],[2.02,0.00]])
A = 3*A     # all loads are multiplied by 3
E_chau = np.sum(A[:,0])  # MWh
q_chau = A[:,0]*1.0e6/hrm    # W
E_clim = np.sum(A[:,1])  # MWh
q_clim = A[:,1]*1.0e6/hrm    # W
q_sol = q_chau+q_clim
qm_ch = max(q_sol)
qm_cl = min(q_sol)
qa = (E_chau + E_clim)*1.0e6/8760   # W
qh_ch = CAP_ch*(COP_ch-1)/COP_ch
qh_cl = -CAP_cl*(COP_cl+1)/COP_cl
#
# time
#
# initial height
H = 110.0
ta = n_years*8760.0
tm = ta + 730.0
tf = tm + nh

#
# debut du design
# Calcul des Tf
#
Tfi_ch = Tfo_ch - qh_ch/(mp*Cpf)
Tfi_cl = Tfo_cl - qh_cl/(mp*Cpf)
Tf_ch = (Tfi_ch+Tfo_ch)/2.0
Tf_cl = (Tfi_cl+Tfo_cl)/2.0
#
# Calul des résistances du sol
#

# calcul de la longeur Chauffage
ok = False
compt = 0
delta = 0.1
zt = np.zeros([nt,2])
qa_cl = min(qa,0)
qa_ch = max(qa,0)
zo = 4
while not ok:
    Bt = b/H
    rr = rb/H
    zot = zo/H
    ib = 0
    for i in range(0,nx):
        for j in range(0,ny):
            zt[ib] = [i*Bt,j*Bt]
            ib = ib+1
    tc = H**2/(9*alhr)
    Fof = tf/tc
    Fo1 = (tf-ta)/tc
    Fo2 = (tf-tm)/tc
#    print Fof,Fo1,Fo2
    gf = compute_g_function(zt,Fof,rbb = rr,zob = zot)
    g1 = compute_g_function(zt,Fo1,rbb = rr,zob = zot)
    g2 = compute_g_function(zt,Fo2,rbb = rr,zob = zot)
#    print gf,g1,g2
    Rsa = (gf-g1)/(2*pi*ks)
    Rsm = (g1-g2)/(2*pi*ks)
    Rsh = (g2)/(2*pi*ks)
#    print Rsa,Rsm,Rsh
    Lch1 = (qa_ch*Rsa)/(To - Tf_ch )
    Lch2 = (qm_ch*Rsm)/(To - Tf_ch )
    Lch3 = (qh_ch*Rsh)/(To - Tf_ch)
    Lch4 = (qh_ch*Rb)/(To - Tf_ch )
    L_ch = Lch1 + Lch2 + Lch3 + Lch4
    # calcul de la longeur Climatisation
    Lcl1 = (qa_cl*Rsa)/(To - Tf_cl )
    Lcl2 = (qm_cl*Rsm)/(To - Tf_cl)
    Lcl3 = (qh_cl*Rsh)/(To - Tf_cl)
    Lcl4 = (qh_cl*Rb)/(To - Tf_cl )
    L_cl = Lcl1 + Lcl2 + Lcl3 + Lcl4
    L = max(L_ch,L_cl)
    print ('The total length is '+ str(L) + ' m')
    Hn = L/nt
    print ('length per borehole is '+ str(Hn) + ' m')
    if abs(Hn-H) < delta:
        ok = True
    else:
        H = Hn
        compt = compt + 1








