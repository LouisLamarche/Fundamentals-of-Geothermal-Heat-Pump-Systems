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

def h_pompe_90r(q,xx):   # q est en l/s, h est en m
    h =   11.725*xx**2  + 0.187*q*xx  -0.085*q**2 -0.00127*q**3/xx
    return h
def rend_pompe_90r(q,xx):   # q est en l/s, h est en m
    h = 26.5  + 16.7*q/xx  -1.38*(q/xx)**2 + 0.004*(q/xx)**3
    return h

p = np.array([-0.0061,-0.179,0.198,20.34])
ald = 0.09 # m2/jr
alhr = ald/24.0
als = alhr/3600
ks = 2.5
To = 10.0


COP_ch = 3.0
# well
# fluid
rhof = 1000
Cpf = 4200.0
muf = 0.0016
nuf = muf/rhof
#
#
# filed configuration
# resistance

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
deb1 = m3s_gpm(gpm1)
deb2 = m3s_gpm(gpm2)
deb3 = m3s_gpm(gpm3)
debv = np.array([deb1,deb2,deb3])
qv = np.array([q1_ch,q2_ch,q3_ch])
qt = q1_ch + q2_ch
DTmax = 3
mmin = qt/(Cpf*DTmax)
Tin = 0
npac = 3
for i in range(0,npac):
    mp1 = debv[i]*rhof
    CCf = mp1*Cpf
    Tout = Tin - qv[i]/CCf
    print(Tin,Tout)
    Tin = (Tout*mp1 + Tin*(mmin-mp1))/mmin
mp = mmin
CCf = mp*Cpf
L1 = 350        # m
Deltep = 300    # Pa/m
SDR = 11
epsilon = 0
g = 9.81
nboucles = 4   # = npac
nbouclesi = 0     # nombre de boucles ou le débit est imposé    = Nin + Nout -1
npipes = 9
npac = 3
VFD = 0
# initalisation des vecteurs et des cellules
qbi = np.zeros(nboucles)
qpipe = np.zeros(npipes)
eppipe = np.zeros(npipes)
Lepipe = np.zeros(npipes)
Kpipe = np.zeros(npipes)
Dpipe = np.zeros(npipes)
# topologie du réseau
mboucle = [np.array([4,-7]),np.array([5,-8]),np.array([6,-9]),np.array([1,2,3,7,8,9])]
qbi[3] = m3s_gpm(gpm)
qbi[0] = m3s_gpm(gpm1)
qbi[1] = m3s_gpm(gpm2)
qbi[2] = m3s_gpm(gpm3)
tube = [3,3,3,1,1,1,3,3,3]
Lpipe = np.array([30,30,350,0,0,0,3,3,3])
Kpipe = np.array([0,0,74.5,2.12e10,2.69,3.68,0,0,0])
SDR = 11
for i in range(0,npipes):
    Dpipe[i] =  sdr_pipe(tube[i],SDR)[0]
Pompe =  [np.array([])] * npipes
#Circ = np.array([-6.59447711e+08, -3.82705350e+06, -3.82977535e+04,  9.47010593e+01])
#Circ = np.array([-9.05132874e+07, -1.18519014e+07, -8.77911946e+03,  2.73200833e+01])
from hpompe_md import pompe_0011,pompe_1935b,pompe_1207
Circ = pompe_0011()
#p = np.array([-1.27008275e+06, -8.50140999e+04,  1.87165130e+02,  1.17264161e+01])
p = np.array([-6.92311603e+07, -1.75599000e+06,  1.67892552e+03,  1.49078296e+02])
#Pompe[2] = 0.8*p*g
Pompe[2] = 1.02*p

Pompe[3] = np.array([-6.24852317e+08, -3.82705350e+06, -4.04181359e+04,  1.05477736e+02])
Pompe[4] = np.array([-7.64915842e+08, -3.82705350e+06, -3.30171824e+04,  7.03863194e+01])
Pompe[5] = np.array([-9.85972791e+08, -3.82705350e+06, -2.56146682e+04,  4.23629023e+01])
if VFD:
    for i in range(0,4):
        Pompe[2][i] = Pompe[2][i]*fac**(i-1)
p =Pompe[0]
my_pipes = []
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nuf)
# calul des débits de conduites
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()
gp = gpm_m3s(qpipes)
print(gp)
exit()

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
hb1 = h_pipe[0] + h_pipe[1]
hb2 = h_pipe[0] + h_pipe[3]
hb3 = h_pipe[2] + h_pipe[3]
print ('Path first heat pump=',hb3)
print ('Path second heat pump= ',hb2)
print ('Path last heat pump= ',hb1)
h_buildinga = max([hb1,hb2,hb3])
print('h building = ' ,h_buildinga)
Dp_field = Deltep*L1
h_field = Dp_field/(rhof*g)
p = np.array([-1.27488209e-03, -8.50140999e-02,  1.86460540e-01,  1.16382932e+01])
px = np.zeros(4)
ql = deb*1000
def fct(x):
    for i in range(0,4):
        px[i] = p[i]*x**(i-1)
    h_pump = np.polyval(px,ql)
    y = h_pump - h_field
    return y
xx = newton(fct,.6)
for i in range(0,4):
    px[i] = p[i]*xx**(i-1)
h_pumpa = np.polyval(px,ql)
h_pumpb = h_pompe_90r(ql,1)
exit()
h_tota = h_buildinga + h_field
h_ft = ft_m(h_tota)

h_pumpa =  h_pompe_90r(ql,1)
W1 = mp*g*h_tota
rend = rend_pompe_90r(ql,1)
W1e = W1/rend*100
print('Fluid power a) = ', W1 ,' W')
print('Electric power a) = ',W1e,' W')
print('Total head  a) = ', h_tota ,' m')
print('Total head  a) = ', h_pumpa + h_buildinga,' m')
h_field_nom = h_field

exit()
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
#
# c
#
h_totc = h_buildinga + h_fieldb
px = np.zeros(4)
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
