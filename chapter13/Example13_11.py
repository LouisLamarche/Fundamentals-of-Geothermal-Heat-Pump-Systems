import numpy as np
from scipy.integrate import nquad,dblquad
from numba import float64,jit
from math import erfc
from matplotlib.pyplot import *
from time import *
import geothermal_md as geo
from CoolProp.CoolProp import *


# main

pi = np.pi
ks = 1.09
als = 7e-7
almin = als*60
alhr = almin*60
D = 0.8
Rs = D/2
Di = 0.024
t = 0.005
Dp = Di + 2*t
kp = 0.39
Rpc = np.log(Dp/Di)/(2*pi*kp)
Trefk = 300
patm = 101000
qv = np.array([10.4,11.3,12.1])
fluid = 'Water'
visc = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
mp = 11*rhof/1000/60
Re = (4*mp)/(pi*Di*visc)           # nombre de Reynold
if (Re>2300):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
        # Dittus-Boetter
else:
    Nud = 3.6
    print('laminar')
hf = Nud*kf/Di
Rpco = 1/(hf*pi*Di)
Rpp = Rpc + Rpco
rp = Dp/2
rpt = rp/Rs
z = 1.5
t =   np.logspace(-.6,2,12)
Fo = alhr*t/Rs**2
nti = len(t)
theta1 = np.zeros(nti)
theta2 = np.zeros(nti)
dd = 2
zt = z/Rs
dtt = dd/Rs
# nombre de rings
ntr = 2

# pitch
p1 = 0.4
Nt1 = int(36/p1)
Nr1 = ntr*Nt1
# longeur tranchée
#Lt = Nt*p + D;print('Lt = ',Lt)
# longeur tuyau
#Lp = Nt*(pi*D+2*p)+D*(pi/2+1);print('Lp = ',Lp)
longu = np.array([500,380,308])
pt1 = p1/Rs
Rb1 = np.zeros([Nr1,2])
for i in range(0,Nt1):
    Rb1[i] = [i*pt1,0]
long1 = (Nt1-1)*pt1
for i in range(Nt1,Nr1):
    Rb1[i] = [long1 - (i-Nt1)*pt1,dtt]
# pitch
p2 = 0.6
Nt2 = int(36/p2)
Nr2 = ntr*Nt2
pt2 = p2/Rs
Rb2 = np.zeros([Nr2,2])
for i in range(0,Nt2):
    Rb2[i] = [i*pt2,0]
long2 = (Nt2-1)*pt2
for i in range(Nt2,Nr2):
    Rb2[i] = [long2 - (i-Nt2)*pt2,dtt]
# pitch
p3 = 0.8
Nt3 = int(36/p3)
Nr3 = ntr*Nt3
pt3 = p3/Rs
Rb3 = np.zeros([Nr3,2])
for i in range(0,Nt3):
    Rb3[i] = [i*pt3,0]
long3 = (Nt3-1)*pt3
for i in range(Nt3,Nr3):
    Rb3[i] = [long3 - (i-Nt3)*pt3,dtt]
theta1 = np.zeros(nti)
theta2 = np.zeros(nti)
theta3 = np.zeros(nti)
r1 = 2.5/Rs
r2 = 50/Rs

tic = time()
for i in range(0,nti):
    g2 = geo.G_function_slinky(Rb1,Fo[i],rpt,zt,pt1,ntr,r1,r2)
    g2b = g2/ks+Rpp*pt1/(2*pi+2*pt1)
    theta1[i] = g2b
toc = time() - tic;print(toc)
tic = time()
for i in range(0,nti):
    g3 = geo.G_function_slinky(Rb2,Fo[i],rpt,zt,pt2,ntr,r1,r2)
    g3b = g3/ks+Rpp*pt2/(2*pi+2*pt2)
    theta2[i] = g3b
toc = time() - tic;print(toc)
for i in range(0,nti):
    g2 = geo.G_function_slinky(Rb3,Fo[i],rpt,zt,pt3,ntr,r1,r2)
    g2 = g2/ks+Rpp*pt3/(2*pi+2*pt3)
    theta3[i] = g2
toc = time() - tic;print(toc)
xx = t/24

data = np.loadtxt('..\\data\\trt1.txt',skiprows = 1)
x1 = data[:,0]
y1 = data[:,1]
data = np.loadtxt('..\\data\\trt2.txt',skiprows = 1)
x2 = data[:,0]
y2 = data[:,1]
data = np.loadtxt('..\\data\\trt3.txt',skiprows = 1)
x3 = data[:,0]
y3 = data[:,1]
p1 = semilogx(xx,theta1,color = 'r',linestyle = 'none',linewidth = 2,marker = 'o',markersize=8,label = 'p = 0.4 m (calculated)')
p2 = semilogx(xx,theta2,color = 'b',linestyle = 'none',linewidth = 2,marker = 'x',markersize=8,label = 'p = 0.6 m (calculated)')
p3 = semilogx(xx,theta3,color = 'g',linestyle = 'none',linewidth = 2,marker = '+',markersize=8,label = 'p = 0.8 m (calculated)')
p4 = semilogx(x1,y1,color = 'r',linestyle = '-',linewidth = 2,label = 'p = 0.4 m (measured)')
p5 = semilogx(x2,y2,color = 'b',linestyle ='-',linewidth = 2,label = 'p = 0.6 m (measured)')
p6 = semilogx(x3,y3,color = 'g',linestyle = '-',linewidth = 2,label = 'p = 0.8 m (measured)')
legend(fontsize = 12)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14

xlabel(' Time ( days)',fontname='Times new Roman',fontsize = fts)
ylabel(r'$\Delta T/q''$',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()


