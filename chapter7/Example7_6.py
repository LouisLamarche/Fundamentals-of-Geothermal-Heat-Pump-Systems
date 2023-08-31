#coding: utf-8
#-------------------------------------------------------------------------------
# Exemple 7.6
#
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse


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
qn1 = newton(fct1,Q)
qn2 = newton(fct2,Q)
h_totn1a = z_tot*qn1**2
qln1 = qn1*1000       # debit en litre/s
h_totn1b = h_pompe(qln1)
print ('Q (keeping the same f) = ' ,qln1, ' l/s')
print ('h total (keeping the same f) = = ' ,h_totn1a, ' m')
h_pi = ft_m(h_totn1a)
gpm = gpm_ls(qln1)
print ('h total (keeping the same f) = ' , h_pi,  ' ft')
print ('gpm (keeping the same f) = ' ,gpm ,' gpm')
h_totn2a = fct_h(qn2)
qln2 = qn2*1000       # debit en litre/s
h_totn2b = h_pompe(qln2)
print ('Q =  (recalculating  f) ' ,qln2, ' l/s')
print ('h total  (recalculating  f) = ' ,h_totn2a, ' m')
print ('h total (recalculating  f) = ' ,  ft_m(h_totn2a),  ' ft')
print ('gpm (recalculating  f) = ' ,gpm_ls(qln2) ,' gpm')


mp1 = qn1*rho
gpmn1 = qn1/m3s_gpm()
P1 = mp1*g*h_totn1a
HP1 = P1/W_hp()
q = np.arange(0.1,2,0.1)         # l/s
h1 =  h_pompe(q)
nn = len(q)
h2 = np.zeros(nn)
for i in range(0,nn):
    h2[i] = fct_h(q[i]/1000)
plot(q,h1,q,h2)
rr = 0.05
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
circle1 = Ellipse((qln1,h_totn1a), rr,rr*ratio, color='r')
ax.add_artist(circle1)
show()
