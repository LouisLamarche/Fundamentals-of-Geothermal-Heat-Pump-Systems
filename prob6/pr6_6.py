from geothermal_md import *
from design_md import *
import  numpy  as np
from  matplotlib.pyplot import *
# sol
ks = 2.95
Cs = 2.8e6
als = ks/Cs
alhr = als*3600
alj = alhr*24
To = 13.0
my_ground = ground(ksoil=ks,alsoil=als,To = To)

# fluide
mp = 4
Cp = 4200.0
CCf = mp*Cp
Tfo_ch =  4.0    # valeur de design
Tfo_cl = 26.0   # valeur de design
#
# puits
#
rb = 0.06
d1 = 0.03
d2 = 0.04
ro = d2/2.0
ri = d1/2.0
xc = rb - ro
kp = 0.6
kc = 1.5
Rcond = np.log(ro/ri)/(2*pi*kp)
h = 1000
Rconv = 1/(pi*d1*h)
Rp = Rcond + Rconv
J = 10
z = np.array([xc,-xc])  # pipes coordinates
Rpcc1,Rpca1 = Rb_linesource(kc,ks,rb,ro,xc)
Rpb1 = Rpcc1 + Rp/2
Rpa1 = Rpca1 + 2*Rp
Rpb,Rpa = Rb_multipole(kc,ks,rb,ro,Rp,J,z)
print(Rpb1,Rpa1)
print(Rpb,Rpa)
nx = 3
ny = 3
nt = nx*ny
b = 6
my_borehole = borehole(nx=nx,ny=ny,rb = rb,dist = b,Rb=Rpb,CCf=CCf,Rint = Rpa)
#
# charges
COP_ch = 4.0
COP_cl = 5.0
CAP_ch = 100000.0
CAP_cl = 90000.0
# kW
E_chau_bat = np.array([21.9 , 18.25, 14.6 , 10.95,  3.65,  0.  ,  0.  ,  0.  ,  0.  ,
       10.95, 14.6 , 21.9 ])
E_clim_bat =  np.array([  0.  ,   0.  ,   0.  ,   0.  ,  -7.3 , -10.95, -14.6 , -14.6 ,
        -7.3 ,   0.  ,   0.  ,   0.  ])
q_chau_bat = E_chau_bat*1e6/730    # W
q_clim_bat = E_clim_bat*1e6/730   # W
q_sol_ch = q_chau_bat*(COP_ch-1)/COP_ch
q_sol_cl = q_clim_bat*(COP_cl+1)/COP_cl
q_sol = q_sol_ch + q_sol_cl
qm_ch = max(q_sol_ch)
qm_cl = min(q_sol_cl)
qa = np.mean(q_sol)  # W
qa2 = sum(q_chau_bat*(COP_ch-1)/COP_ch +q_clim_bat*(COP_cl+1)/COP_cl)/12
qh_ch = CAP_ch*(COP_ch-1)/COP_ch
qh_cl = -CAP_cl*(COP_cl+1)/COP_cl
#
# temps
#
n_years  = 20
nh = 6

my_borehole = borehole(nx=nx,ny=ny,rb = rb,dist = b,Rb=Rpb,CCf=CCf,Rint = Rpa)
param_conception_eed = eed_params(init_month = 0,q_months = q_sol,qh_heating=qh_ch \
        ,qh_cooling=qh_cl,Tfo_heating =Tfo_ch,Tfo_cooling \
        = Tfo_cl,n_years = n_years, n_bloc = nh,flag_inter=True)
second_design = borefield(params = param_conception_eed,ground = my_ground,borehole = my_borehole)
L2 =  second_design.Compute_L_eed()
print('L = ',L2)
print('H = ',L2/nt)
print('Rpb = ',Rpb)
print('Rps = ',second_design.Rbs)
print('Tfi (heating) = ',second_design.Tfi_he)
print('Tfi (cooling) = ',second_design.Tfi_co)
print('Tf (heating) = ',(second_design.Tfi_he+Tfo_ch)/2)
print('Tf (cooling) = ',(second_design.Tfi_co+Tfo_cl)/2)
print('L (heating) = ',second_design.L_he)
print('L (cooling) = ',second_design.L_co)



#
