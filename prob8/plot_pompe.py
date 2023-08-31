from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from copy import *
from conversion_md import *
from scipy.optimize import newton

def h_pompe_90(q,xx):   # q est en l/s, h est en m
    h =   20.34*xx**2  + 0.198*q*xx  -0.179*q**2 -0.0061*q**3/xx
    return h
p = np.array([-0.0061,-0.179,0.198,20.34])
def rend_pompe_90(q,xx):   # q est en l/s, h est en m
    h = 17.96  + 26.28*q/xx  -4.07*(q/xx)**2 + 0.157*(q/xx)**3
    return h



hft = 38.12
gpm = 72
h_totc = m_ft(hft)
qlb = ls_gpm(gpm)
gpv = arange(0,80)
def fct(x):
    h_pump = h_pompe_90(qlb,x)
    y = h_pump - h_totc
    return y
xx = newton(fct,.8)
px = zeros(4)
for i in range(0,4):
    px[i] = p[i]*xx**(i-1)
n = len(gpv)
h = zeros(n)
h2 = zeros(n)
renv = zeros(n)
for i in range(0,n):
    q = ls_gpm(gpv[i])
    hh = h_pompe_90(q,xx)
    hx = polyval(px,q)
    renv[i] = rend_pompe_90(q,xx)
    h[i] = ft_m(hh)
    h2[i] = ft_m(hx)
figure(1)
plot(gpv,h,gpv,h2,'x')
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
rr = 2.5
rrn = 2.5
dx = rrn/2
dy = rrn*ratio/2
circle1 = Ellipse((gpm,hft), rr,rr*ratio, color='r')
ax.add_artist(circle1)
figure(2)
ren =rend_pompe_90(qlb,xx)
plot(gpv,renv)
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
rr = 2.5
rrn = 2.5
dx = rrn/2
dy = rrn*ratio/2
circle1 = Ellipse((gpm,ren), rr,rr*ratio, color='r')
ax.add_artist(circle1)
show()