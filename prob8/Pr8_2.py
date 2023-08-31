#
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
import numpy as np
from  matplotlib.pyplot import *
from scipy.optimize import newton
#
# data
#
# soil

def h_pompe_90(q,xx):   # q est en l/s, h est en m
    h =   15.5*xx**2  + 0.177*q*xx  -0.183*q**2 - 0.0074*q**3/xx
    return h

def rend_pompe_90(q,xx):   # q est en l/s, h est en m
    h = 17.96  + 30.3*q/xx  - 5.45*(q/xx)**2 + 0.243*(q/xx)**3
    return h


# fluid
rhof = 1000
Cpf = 4200.0
muf = 0.0016
nuf = muf/rhof
#
# loads MWh
#
q1= 110 # MBTUh
q2= 90 # MBTUh
q3= 75 # MBTUh
gpm1 = 30  # gpm
gpm2 = 24  # gpm
gpm3 = 18  # gpm
gpm_hpv = np.array([gpm1,gpm2,gpm3])
psi_hpv = np.array([1.6,1.3,1.0])
gpm = gpm1 + gpm2 + gpm3
q1W = 1000*W_BTUhr(q1)
q2W = 1000*W_BTUhr(q2)
q3W = 1000*W_BTUhr(q3)
deb1 = m3s_gpm(gpm1)  #
deb2 = m3s_gpm(gpm2)
deb3 = m3s_gpm(gpm3)
deb = (deb1+deb2+deb3)
mp = (deb1+deb2+deb3)*rhof
mp1 = deb1*rhof
mp2 = deb2*rhof
mp3 = deb3*rhof
CCf = mp*Cpf
L1 = 350        # m
Deltep = 300    # Pa/m
SDR = 11
epsilon = 0
g = 9.81
#
# a)
#
h_hpv = mw_psi(psi_hpv)
tube = np.array([3,3,3,3,3])
Lpipe  = np.array([3,30,3,30,3])
gpmv = np.array([gpm - gpm1,gpm,gpm-gpm2,gpm,gpm-gpm3])
h_pipe = np.ones(5)
h_building = 0
for i in range(0,5):
    di,do = sdr_pipe(tube[i],SDR)
    A1 = pi*di**2/4.0
    debiti = m3s_gpm(gpmv[i])
    u1 = debiti/A1
    Re = di*u1/nuf
    ed = epsilon/di
    f1 = Colebrook(Re,ed)
    Le = 0
    h_pipe[i] = f1*((Lpipe[i])/di)*u1**2/(2*g) # Dp en metres
    print('pipe number #',i+1,' L equivalent = ','%.2f' % Le,' h total = ',h_pipe[i])
    h_building = h_building + h_pipe[i]
Dp_field = Deltep*L1
h_field = Dp_field/(rhof*g)
h_tota = h_building + h_field
h_ft = ft_m(h_tota)
ql = deb*1000
h_pumpa =  h_pompe_90(ql,1)
W_pumpa = mp*g*h_tota
h_tot_hp1 = h_pipe[0] + h_hpv[0]
h_tot_hp2 = h_pipe[2] + h_hpv[1]
h_tot_hp3 = h_pipe[4] + h_hpv[2]
W_pump_hp1 = mp1*g*h_tot_hp1
rend1 = 0.321*W_pump_hp1**0.115
W_pump_hp2 = mp2*g*h_tot_hp2
rend2 = 0.321*W_pump_hp2**0.115
W_pump_hp3 = mp3*g*h_tot_hp3
rend3 = 0.321*W_pump_hp3**0.115
W_fluid = W_pumpa +W_pump_hp1 + W_pump_hp2 + W_pump_hp3
rend = rend_pompe_90(ql,1)
rend = rend/100
print('Fluid power filed = ',  W_pumpa ,' W')
print('Fluid power HP1 = ',  W_pump_hp1 ,' W')
print('Fluid power HP1 = ',  W_pump_hp2 ,' W')
print('Fluid power HP1 = ',  W_pump_hp3 ,' W')
W_elec = W_pumpa/rend + W_pump_hp1/rend1 +W_pump_hp2/rend2 +W_pump_hp3/rend3
print('Fluid power a) = ', W_fluid ,' W')
print('Electric power a) = ',W_elec,' W')
print('Total head  a) = ', h_tota ,' m')
print('Total head  a) = ', h_pumpa ,' m')
h_field_nom = h_field

# b

Tin1 = 0
Tout1 = Tin1 - q1W/(mp1*Cpf)
Tin2 = (Tout1*mp1 + Tin1*(mp-mp1))/mp
print('Tout 1 = ',Tout1)
print('Tin 2 = ',Tin2)
Tout2 = Tin2 - q2W/(mp2*Cpf)
Tin3 = (Tout2*mp2 + Tin2*(mp-mp2))/mp
print('Tout = ',Tout2)
print('Tin = ',Tin3)
Tout3 = Tin3 - q3W/(mp3*Cpf)
Tin4 = (Tout3*mp3 + Tin3*(mp-mp3))/mp
print('Tout = ',Tout3)
print('Tin = ',Tin4)
Tin5 = Tin1 - (q1W + q2W + q3W)/(mp*Cpf)
print('Tin field = ',Tin5)


# c
def fctb(debn):
    debv = np.array([debn,debn,debn - deb2,debn,debn-deb3])
    h_fieldb = h_field_nom*(debn/deb)**2
    h_buildingb = 0
    for i in range(0,5):
        di,do = sdr_pipe(tube[i],SDR)
        A1 = pi*di**2/4.0
        debiti = debv[i]
        u1 = debiti/A1
        Re = di*u1/nuf
        ed = epsilon/di
        f1 = Colebrook(Re,ed)
        h_p = f1*((Lpipe[i])/di)*u1**2/(2*g) # Dp en metres
        h_buildingb = h_buildingb + h_p
    h_totb = h_buildingb + h_fieldb
    qlb = debn*1000
    h_pumpb =  h_pompe_90(qlb,1)
    return h_totb - h_pumpb
debb = newton(fctb,deb)
h_fieldb = h_field_nom*(debb/deb)**2
h_pumpb =  h_pompe_90(debb*1000,1)
W_pumpb = debb*rhof*g*h_pumpb
gpmb = gpm_m3s(debb)
print('gpm b) = ',gpmb)
W_pump_hp2 = mp2*g*h_tot_hp2
rend2 = 0.321*W_pump_hp2**0.115
W_pump_hp3 = mp3*g*h_tot_hp3
rend3 = 0.321*W_pump_hp3**0.115
W_fluid = W_pumpb  + W_pump_hp2 + W_pump_hp3
rend = rend_pompe_90(debb*1000,1)
rend = rend/100
print('Fluid power filed = ',  W_pumpb ,' W')
W_elec = W_pumpb/rend  +W_pump_hp2/rend2 +W_pump_hp3/rend3
print('Fluid power b) = ', W_fluid ,' W')
print('Electric power b) = ',W_elec,' W')
