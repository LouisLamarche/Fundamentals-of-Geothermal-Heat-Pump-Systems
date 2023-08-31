#-------------------------------------------------------------------------------
# Exemple 7.7
#
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse

pi = np.pi
g = 9.81
epsilon = 0
fluide1 = 'water'
fluide2 = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 20 % volume
fluide3 = 'INCOMP::APG-40%'  # ASHRAE propylene glycl 40 % volume
cas = 1
if cas ==1:
    fluide = fluide1
elif cas ==2:
    fluide = fluide2
else:
    fluide = fluide3
patm = 101.325*1000.0
#
Cv = 1.08  # l/s/kPa**1/2
K = 2.5
L = 100.0
D = 0.04
Le = 2.0
T_refk = 278
Ql = 1.5       # debit el l/s
Q = Ql/1000.0   # debit en m3/s
A = pi*D**2/4.0
u = Q/A
ed  = epsilon/D
rho1 =PropsSI('D','T',T_refk,'P',patm,'water')
mu1 = PropsSI('viscosity','T',T_refk,'P',patm,'water')
nu1 = mu1/rho1
Re1 = D*u/nu1
f1 = Colebrook(Re1,ed)
rhow =PropsSI('D','T',16+273.15,'P',patm,'water')
rho =PropsSI('D','T',T_refk,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk,'P',patm,fluide)
Sg = rho/rhow
nu = mu/rho
Re = D*u/nu
f = Colebrook(Re,ed)
corr = f/f1
z_pipe = f*(L/D)/(A**2*2*g) # Dp en metres
z_sing1 = K/(A**2*2*g)*corr
z_sing2 = f*(Le/D)/(A**2*2*g)
DPPa = 1000*Sg*(1e3/Cv)**2
z_valve = DPPa/(rho*g)*corr # m de fluide
z_sing = z_sing1 + z_sing2
z_tot = z_pipe + z_sing + z_valve
p = np.array([- 0.01588,-0.15556,-0.04935,5.4952])
def h_pompe(q):   # q est en l/s
    h = np.polyval(p,q)                            # h est en metres
    return h
def fct1(q):
    hs = z_tot*q**2
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp = h_pompe(ql)
    y = hs - hp
    return y
def fct_h(q):
    u = q/A
    Re1 = D*u/nu1
    f1 = Colebrook(Re1,ed)
    Re = D*u/nu
    f = Colebrook(Re,ed)
    corr = f/f1
    z_pipe = f*(L/D)/(A**2*2*g) # Dp en metres
    z_sing1 = K/(A**2*2*g)*corr
    z_sing2 = f*(Le/D)/(A**2*2*g)
    z_valve = DPPa/(rho*g)*corr # m de fluide
    z_sing = z_sing1 + z_sing2
    z_tot = z_pipe + z_sing + z_valve
    hs = z_tot*q**2
    return hs

def fct2(q):
    hs = fct_h(q)
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp = h_pompe(ql)
    y = hs - hp
    return y
qn = newton(fct2,Q)
qln = qn*1000       # debit en litre/s
h_totn = h_pompe(qln)
#
# 80 % reduction
#
fac = 0.8
Dpompe = m_ft(4.25/12)
om = 1750*2*pi/60
om2 = fac*om
Cq1 = qn/(om*Dpompe**3)
Ch1 = g*h_totn/(om**2*Dpompe**2)
qr1 = Cq1*om2*Dpompe**3
hr1 = Ch1*(om2**2*Dpompe**2)/g
gpmr = gpm_m3s(qr1)
hr1_pi = ft_m(hr1)
print ('Q (80 % reduction ) = ' ,qr1*1000, ' l/s')
print ('gpm (80 % reduction )= ' ,gpmr, ' gpm')
print ('h total  (80 % reduction )=' ,hr1, ' m')
print ('h total (80 % reduction ) = ', hr1_pi,' ft')
px = np.zeros(4)
for i in range(0,4):
    px[i] = p[i]*fac**(i-1)
def fct1n(q):
    hs = z_tot*q**2
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp = np.polyval(px,ql)
    y = hs - hp
    return y

def fct2n(q):
    hs = fct_h(q)
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp =  np.polyval(px,ql)
    y = hs - hp
    return y
qr2 = newton(fct1n,Q)
hr2a = z_tot*qr2**2
hr2b = np.polyval(px,qr2*1000)
gpmr2 = gpm_m3s(qr2)
hr2_pi = ft_m(hr2a)
print ('Q (80 % reduction ) = ' ,qr2*1000, ' l/s')
print ('gpm (80 % reduction )= ' ,gpmr2, ' gpm')
print ('h total  (80 % reduction )=' ,hr2a, ' m')
print ('h total  (80 % reduction )=' ,hr2b, ' m')
print ('h total (80 % reduction ) = ', hr2_pi,' ft')
print('')
print('Taking into account variation of f')
print('')
qr3 = newton(fct2n,Q)
hr3a = fct_h(qr3)
hr3b = np.polyval(px,qr3*1000)
gpmr3 = gpm_m3s(qr3)
hr3_pi = ft_m(hr3a)
print ('Q (80 % reduction ) = ' ,qr3*1000, ' l/s')
print ('gpm (80 % reduction )= ' ,gpmr3, ' gpm')
print ('h total  (80 % reduction )=' ,hr3a, ' m')
print ('h total  (80 % reduction )=' ,hr3b, ' m')
print ('h total (80 % reduction ) = ', hr3_pi,' ft')

q = np.arange(0.1,1.8,0.1)
h1 =  h_pompe(q)
h2 = z_tot*(q/1000)**2
h3 =  np.polyval(px,q)
nn = len(q)
h4 = np.zeros(nn)
for i in range(0,nn):
    h4[i] = fct_h(q[i]/1000)
plot(q,h1,q,h2,q,h3,q,h4)
rr = 0.05
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
circle1 = Ellipse((qln,h_totn), rr,rr*ratio, color='r')
ax.add_artist(circle1)
circle1 = Ellipse((qr2*1000,hr2a), rr,rr*ratio, color='r')
ax.add_artist(circle1)
circle2 = Ellipse((qr3*1000,hr3a), rr,rr*ratio, color='r')
ax.add_artist(circle2)
xticks(fontsize=14)
yticks(fontsize=14)
show()
