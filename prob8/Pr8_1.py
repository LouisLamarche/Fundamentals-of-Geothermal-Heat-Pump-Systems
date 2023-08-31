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
    h =   20.34*xx**2  + 0.198*q*xx  -0.179*q**2 -0.0061*q**3/xx
    return h

def rend_pompe_90(q,xx):   # q est en l/s, h est en m
    h = 17.96  + 26.28*q/xx  -4.07*(q/xx)**2 + 0.157*(q/xx)**3
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
deb1 = m3s_gpm(gpm1)
deb2 = m3s_gpm(gpm2)
deb3 = m3s_gpm(gpm3)
deb = (deb1+deb2+deb3)
mp = (deb1+deb2+deb3)*rhof
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
gpmv = np.array([gpm - gpm1,gpm - gpm1 - gpm2,gpm1,gpm1+gpm2])
tube = np.array([2,1.5,1.5,2])
Lpipe = 30*np.ones(4)
h_pipe = np.ones(4)
for i in range(0,4):
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
hb1 = h_pipe[0] + h_pipe[1] + h_hpv[2] + m_ft(3.5)
hb2 = h_pipe[0] + h_pipe[3] + h_hpv[1] + m_ft(3.5)
hb3 = h_pipe[2] + h_pipe[3] + h_hpv[0] + m_ft(3.5)


print ('Path first heat pump=',hb3)
print ('Path second heat pump= ',hb2)
print ('Path last heat pump= ',hb1)
h_buildinga = max([hb1,hb2,hb3])
print('h building = ' ,h_buildinga)
Dp_field = Deltep*L1
h_field = Dp_field/(rhof*g)
h_tota = h_buildinga + h_field
h_ft = ft_m(h_tota)
ql = deb*1000
h_pumpa =  h_pompe_90(ql,1)
W1 = mp*g*h_tota
W11 = mp*g*h_buildinga
W13 = mp*g*h_field
rend = rend_pompe_90(ql,1)
W1e = W1/rend*100
print('Fluid power a) = ', W1 ,' W')
print('Fluid power build a) = ', W11 ,' W')
print('Fluid power field a) = ', W13 ,' W')
print('Electric power a) = ',W1e,' W')
print('Total head  a) = ', h_tota ,' m')
print('Total head  a) = ', h_pumpa ,' m')
h_field_nom = h_field


#
# b
#
gpmb = gpm2 + gpm3
debb = m3s_gpm(gpmb)
mpb = debb*rhof
gpmvb = np.array([gpmb,gpmb - gpm2 ,0,gpm2])
h_pipeb = np.ones(4)
for i in range(0,4):
    di,do = sdr_pipe(tube[i],SDR)
    A1 = pi*di**2/4.0
    debiti = m3s_gpm(gpmvb[i])
    u1 = debiti/A1
    Re = di*u1/nuf
    ed = epsilon/di
    f1 = Colebrook(Re,ed)
    Le = 0
    h_pipeb[i] = f1*((Lpipe[i])/di)*u1**2/(2*g) # Dp en metres
    print('pipe number #',i+1,' L equivalent = ','%.2f' % Le,' h total = ',h_pipe[i])
h_fieldb = h_field_nom*(gpmb/gpm)**2
def fctb(h_ctrl):
    hb1 = h_pipe[0] + h_pipe[1] + h_hpv[2] + h_ctrl
    hb2 = h_pipe[0] + h_pipe[3] + h_hpv[1] +  + h_ctrl
    h_buildingb = max([hb1,hb2])
    h_totb = h_buildingb + h_fieldb
    qlb = debb*1000
    h_pumpb =  h_pompe_90(qlb,1)
    return h_totb - h_pumpb
h_ctr = newton(fctb,1)
hb1 = h_pipe[0] + h_pipe[1] + h_hpv[2] + h_ctr
hb2 = h_pipe[0] + h_pipe[3] + h_hpv[1] +   h_ctr
h_buildingb = max([hb1,hb2])
h_totb = h_buildingb + h_fieldb
qlb = debb*1000
h_pumpb =  h_pompe_90(qlb,1)
W1b = mpb*g*h_totb
rendb = rend_pompe_90(qlb,1)
W1eb = W1b/rendb*100
print('Fluid power b) = ', W1b ,' W')
print('Electric power b) = ',W1eb,' W')
print('Total head  b) = ', h_totb ,' m')
print('Total head  b) = ', h_pumpb ,' m')
C
#
# c
#
h_totc = h_buildinga + h_fieldb
def fct(x):
    h_pump = h_pompe_90(qlb,x)
    y = h_pump - h_totc
    return y
xx = newton(fct,.6)
h_pumpc =  h_pompe_90(qlb,xx)
W1c = mpb*g*h_totc
rendc = rend_pompe_90(qlb,xx)
W1ec = W1c/rendc*100
print('Fluid power c) = ', W1c ,' W')
print('Electric power c) = ',W1ec,' W')
print('Total head  c) = ', h_totc ,' m')
print('Total head  c) = ', h_pumpc ,' m')
