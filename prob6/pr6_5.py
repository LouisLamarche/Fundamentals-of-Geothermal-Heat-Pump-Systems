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
#xc = 0.035
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
Hi = 175
mp1 = mp/9
eta = Hi/(mp1*Cp*np.sqrt(Rpb*Rpa))
Rbs1 = Rpb*eta/np.tanh(eta)
Rbs2 = Rpb*(1 + eta**2/3)
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
A = np.array([[30,0.00],[25, 0.00],[20,0.0],[15,0],[5,-10],[0,-15], \
            [0.0,-20],[0.00,-20],[0.0,-10],[15,0.0],[20,0.0],[30,0.0]])
q_chau_bat = A[:,0]*1000    # W
q_clim_bat = A[:,1]*1000    # W
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
param_conception_ashrae = ashrae_params(qa=qa,qm_heating = qm_ch, qm_cooling=qm_cl,\
    qh_heating=qh_ch,qh_cooling=qh_cl,Tfo_heating =Tfo_ch,\
    Tfo_cooling = Tfo_cl,n_years = n_years, n_bloc = nh,flag_inter=True,flag_Tp = 'FLS')
first_design = borefield(params = param_conception_ashrae,ground = my_ground,borehole = my_borehole)
L1 =  first_design.Compute_L_ashrae()
print('L = ',L1)
print('H = ',L1/nt)
print('Tp = ',first_design.Tp)
print('Rpb = ',Rpb)
print('Rps = ',first_design.Rbs)
print('Tfi (heating) = ',first_design.Tfi_he)
print('Tfi (cooling) = ',first_design.Tfi_co)
print('Tf (heating) = ',(first_design.Tfi_he+Tfo_ch)/2)
print('Tf (cooling) = ',(first_design.Tfi_co+Tfo_cl)/2)
print('L (heating) = ',first_design.L_he)
print('L (cooling) = ',first_design.L_co)
