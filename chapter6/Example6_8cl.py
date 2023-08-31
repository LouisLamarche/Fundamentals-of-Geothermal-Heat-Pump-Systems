
# Exemple 6.8 written by Louis lamarche 3 october 2017
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from design_md import *
#
# data
#
# soil
ald = 0.09 # m2/jr
alhr = ald/24.0
als = alhr/3600
ks = 2.5
To = 10.0

my_ground = ground(ksoil=ks,alsoil=als,To = To)
#
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
# filed configuration
nx = 3
ny = 2
b = 6.0
b = 6.0
nt = nx*ny
mp = mp1*nt
zot = 0.04
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
CCf = mp*Cpf
my_borehole = borehole(nx=nx,ny=ny,rb = rb,dist = b,Rb=Rb,CCf=CCf)
param_conception_ashraea = ashrae_params(qa=qa,qm_heating = qm_ch, qm_cooling=qm_cl,\
    qh_heating=qh_ch,qh_cooling=qh_cl,Tfo_heating =Tfo_ch,\
    Tfo_cooling = Tfo_cl,n_years = n_years, n_bloc = nh,flag_Tp = 'ASHRAE')
first_design = borefield(params = param_conception_ashraea,ground = my_ground,borehole = my_borehole)
L1 =  first_design.Compute_L_ashrae()
print('L = ',L1)
print('Tp = ',first_design.Tp)
param_conception_ashraeb = ashrae_params(qa=qa,qm_heating = qm_ch, qm_cooling=qm_cl,\
    qh_heating=qh_ch,qh_cooling=qh_cl,Tfo_heating =Tfo_ch,\
    Tfo_cooling = Tfo_cl,n_years = n_years, n_bloc = nh,flag_Tp = 'FLS')
second_design = borefield(params = param_conception_ashraeb,ground = my_ground,borehole = my_borehole)
L2 =  second_design.Compute_L_ashrae()
print('L = ',L2)
print('Tp = ',second_design.Tp)



