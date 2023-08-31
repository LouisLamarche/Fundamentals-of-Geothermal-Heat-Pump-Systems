
# Exemple 6.10 written by Louis lamarche 3 october 2017
#
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
import numpy as np
from  matplotlib.pyplot import *
from design_md import *
from scipy.optimize import newton
#
# data
#
# soil

def h_pompe_90(q):   # q est en l/s, h est en m
    h =   20.34  + 0.198*q  -0.179*q**2 -0.0061*q**3
    return h

def rend_pompe_90(q):   # q est en l/s, h est en m
    h = 17.96  + 26.28*q  -4.07*q**2 + 0.157*q**3
    return h
p = np.array([-0.0061,-0.179,0.198,20.34])
ald = 0.09 # m2/jr
alhr = ald/24.0
als = alhr/3600
ks = 2.5
To = 10.0

my_ground = ground(ksoil=ks,alsoil=als,To = To)

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
rhof = 1000
Cpf = 4200.0
muf = 0.0016
nuf = muf/rhof
#
# design data
Tfo_ch = 0.0    # valeur de design
Tfo_cl = 25.0   # valeur de design
n_years = 10
nh = 4          # hourly block
#
# filed configuration
# resistance

Rb = 0.1
#
# loads MWh
#
q1= 110 # MBTUh
q2= 90 # MBTUh
q3= 75 # MBTUh
gpm1 = 30  # gpm
gpm2 = 24  # gpm
gpm3 = 18  # gpm
psi_hpv = np.array([1.6,1.3,1.0])
gpm = gpm1 + gpm2 + gpm3
q1W = 1000*W_BTUhr(q1)
q2W = 1000*W_BTUhr(q2)
q3W = 1000*W_BTUhr(q3)
q1_ch = q1W*(COP_ch-1)/COP_ch
q2_ch = q2W*(COP_ch-1)/COP_ch
q3_ch = q3W*(COP_ch-1)/COP_ch
q_ch = (q1_ch+q2_ch+q3_ch)
q1_cl = -q1W*(COP_cl+1)/COP_cl
q2_cl = -q2W*(COP_cl+1)/COP_cl
q3_cl = -q3W*(COP_cl+1)/COP_cl
q_cl = (q1_cl+q2_cl+q3_cl)
deb1 = m3s_gpm(gpm1)
deb2 = m3s_gpm(gpm2)
deb3 = m3s_gpm(gpm3)
deb = (deb1+deb2+deb3)
mp = (deb1+deb2+deb3)*rhof
CCf = mp*Cpf
Fo = alhr*6/rb*2
TLSch = q_ch*G_function(Fo)/ks
TLScl = q_cl*G_function(Fo)/ks
Tfi_ch = Tfo_ch+q_ch/CCf
Tfi_cl = Tfo_cl+q_cl/CCf
Tf_ch = (Tfi_ch + Tfo_ch)/2
Tf_cl = (Tfi_cl + Tfo_cl)/2
La =  (TLSch + q_ch*Rb)/(To - Tf_ch)
Lb =  (TLScl + q_cl*Rb)/(To - Tf_cl)
print('La = ',La,Lb)
L1 = 350
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
h_building = max([hb1,hb2,hb3])
print('h building = ' ,h_building)
Dp_field = 300*L1
h_field = Dp_field/(rhof*g)
h_tot = h_building + h_field
h_ft = ft_m(h_tot)
ql = deb*1000
h_pump =  h_pompe_90(ql)
W1 = mp*g*h_tot
rend = rend_pompe_90(ql)
W1e = W1/rend*100
print(W1,W1e)
print(h_tot,h_pump)


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
    h_pumpb =  h_pompe_90(qlb)
    return h_totb - h_pumpb
h_ctr = newton(fctb,1)
hb1 = h_pipe[0] + h_pipe[1] + h_hpv[2] + h_ctr
hb2 = h_pipe[0] + h_pipe[3] + h_hpv[1] +   h_ctr
h_buildingb = max([hb1,hb2])
h_totb = h_buildingb + h_fieldb
qlb = debb*1000
h_pumpb =  h_pompe_90(qlb)
h_pumpbb =  np.polyval(p,qlb)
W1b = mpb*g*h_totb
rendb = rend_pompe_90(qlb)
W1eb = W1b/rendb*100
print('b',W1b,W1eb)
print(h_totb,h_pumpb,h_pumpbb)
#
# c
#
h_totc = h_building + h_fieldb
px = np.zeros(4)
def fct(x):
    for i in range(0,4):
        px[i] = p[i]*x**(i-1)
    h_pump = np.polyval(px,qlb)
    y = h_pump - h_totc
    return y
xx = newton(fct,.6)
for i in range(0,4):
    px[i] = p[i]*xx**(i-1)
h_pumpc =  np.polyval(px,qlb)
W1c = mpb*g*h_totc
rendc = rend_pompe_90(qlb/xx)
W1ec = W1c/rendc*100
print('c',W1c,W1ec)
print(h_totc,h_pumpc)

