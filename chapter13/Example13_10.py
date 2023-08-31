from numpy import *
from scipy.integrate import nquad,dblquad
from numba import float64,jit
from math import erfc
from matplotlib.pyplot import *
from time import *
import geothermal_md as geo


# main


cas = 'b'
if cas == 'a':
    ylab = r'$g_{slinky}$'
else:
    ylab = r'$g_{slinky}\tilde{p}/(2\pi+\tilde{p})$'
ks = 1.12
als = 4.33e-7
almin = als*60
alhr = almin*60
D = 0.8
Rs = D/2
Di = 0.024
t = 0.005
Dp = Di + 2*t
kp = 0.4
Rpc = log(Dp/Di)/(2*pi*kp)
Rpco = 1/(1000*pi*Di)
Rpp = Rpc + Rpco
rp = Dp/2
rpt = rp/Rs
z = 1.5
t =   logspace(-2,5,20)
Fo = alhr*t/Rs**2
nti = len(t)
theta1 = zeros(nti)
theta2 = zeros(nti)
dd = 2
zt = z/Rs
dtt = dd/Rs
# nombre de trench
ntr = 2
# pitch
p1 = 0.4
Nt1 = int(36/p1)
Nr1 = ntr*Nt1
pt1 = p1/Rs
Rb1 = zeros([Nr1,2])
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
Rb2 = zeros([Nr2,2])
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
Rb3 = zeros([Nr3,2])
for i in range(0,Nt3):
    Rb3[i] = [i*pt3,0]
long3 = (Nt3-1)*pt3
for i in range(Nt3,Nr3):
    Rb3[i] = [long3 - (i-Nt3)*pt3,dtt]
theta1 = zeros(nti)
theta2 = zeros(nti)
theta3 = zeros(nti)
r1 = 2.5/Rs  # if r > r1 g = g(approx)
r2 = 120/Rs  # if r > r2 g = 0
tic = time()
for i in range(0,nti):
    g2 = geo.calcul_fonction_g_xiong(Rb1,Fo[i],rpt,zt,Nt1,ntr,r1,r2)
    if cas == 'b':
        g2 = g2*pt1/(2*pi+2*pt1)
    theta1[i] = g2
toc = time() - tic;print(toc)
tic = time()
for i in range(0,nti):
    g2 = geo.calcul_fonction_g_xiong(Rb2,Fo[i],rpt,zt,Nt2,ntr,r1,r2)
    if cas == 'b':
        g2 = g2*pt2/(2*pi+2*pt2)
    theta2[i] = g2
toc = time() - tic;print(toc)
tic = time()
for i in range(0,nti):
    g2 = geo.calcul_fonction_g_xiong(Rb3,Fo[i],rpt,zt,Nt3,ntr,r1,r2)
    if cas == 'b':
        g2 = g2*pt3/(2*pi+2*pt3)
    theta3[i] = g2
toc = time() - tic;print(toc)
xx = log10(t)

p1 = plot(xx,theta1,color = 'k',linewidth = 2,marker = 'o',markersize=8,label = 'p = 0.4 m ')
p2 = plot(xx,theta2,color = 'k',linewidth = 2,marker = 'x',markersize=8,label = 'p = 0.6 m')
p3 = plot(xx,theta3,color = 'k',linewidth = 2,marker = '+',markersize=8,label = 'p = 0.8 m')
legend(fontsize = 12)
#title('$r_p/R = 0.25,h/R = 0.5$',fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14

xlabel(' $log_{10}(t ( hr))$',fontname='Times new Roman',fontsize = fts)
ylabel(ylab,fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()

