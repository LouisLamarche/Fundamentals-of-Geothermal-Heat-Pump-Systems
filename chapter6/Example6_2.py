
# Exemple 6.2 written by Louis lamarche
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
#
# data
#
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])              # hours in each month
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])   # cumulate hours

cas = 'a'
# soil
als = 0.09 # m2/jr
alhr = als/24.0
ks = 2.5
To = 10.0
COP_ch = 3.0
COP_cl = 5.0
modi = 0  # Ex 6.2 modified modi =1
if modi ==1:
    COP_cl = 3
CAP_ch = 12000.0
CAP_cl = 12000.0
# well
rb = 0.15/2.0
Rpb = 0.13
kg = 1.7
do = 0.033
di = 0.027
ro = do/2.0
ri = di/2.0
xc = rb/3
kp = 0.4
# fluid
rhof = 1030
mp = 0.6/1000.0*rhof
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
nh = 4          # hourly block
# resistance
Re = 4*mp/(pi*di*muf)
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
Rb2 = Rb_Sharqawi(kg,rb,ro,xc) + Rp/2.0
Rb3,Ra3 = Rb_linesource(kg,ks,rb,ro,xc) + Rp/2.0
if cas == 'a':
    Rb = Rb1
elif cas == 'b':
    Rb = Rb2
else:
    Rb = Rb3
#
# loads MWh
#
hrm = 730
A = np.array([[2.17,0.00],[2.07, 0.00],[1.750,0.000],[1.39,0.00],[0.9,-1.2],[0.00,-1.6],[0.00,-2.00],[0.00,-1.8],[0.85,-0.80],[1.22,-0.40],[1.64,0.00],[2.02,0.00]])
E_chau = sum(A[:,0])  # MWh
q_chau = A[:,0]*1.0e6/hrm    # W
E_clim = sum(A[:,1])  # MWh
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
ta = n_years*8760.0
tm = ta + 730.0
tf = tm + nh
Fof = alhr*tf/rb**2
Fo1 = alhr*(tf-ta)/rb**2
Fo2 = alhr*(tf-tm)/rb**2
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
Gf = G_function(Fof)
G1 = G_function(Fo1)
G2 = G_function(Fo2)
Rsa = (Gf-G1)/ks
Rsm = (G1-G2)/ks
Rsh = (G2)/ks
# calcul de la longeur Chauffage
Tp = 0.0
Lch1 = (qa*Rsa)/(To - Tf_ch - Tp)
Lch1 = max(Lch1,0)
Lch2 = (qm_ch*Rsm)/(To - Tf_ch - Tp)
Lch3 = (qh_ch*Rsh)/(To - Tf_ch- Tp)
Lch4 = (qh_ch*Rb)/(To - Tf_ch - Tp)
L_ch = Lch1 + Lch2 + Lch3 + Lch4
# calcul de la longeur Climatisation
Lcl1 = (qa*Rsa)/(To - Tf_cl - Tp)
Lcl1 = max(Lcl1,0)
Lcl2 = (qm_cl*Rsm)/(To - Tf_cl - Tp)
Lcl3 = (qh_cl*Rsh)/(To - Tf_cl- Tp)
Lcl4 = (qh_cl*Rb)/(To - Tf_cl - Tp)
L_cl = Lcl1 + Lcl2 + Lcl3 + Lcl4
L = max(L_ch,L_cl)
print ('The heating length is '+ str(L_ch) + ' m')
print ('The cooling length is '+ str(L_cl) + ' m')
print ('The total length is '+ str(L) + ' m')









