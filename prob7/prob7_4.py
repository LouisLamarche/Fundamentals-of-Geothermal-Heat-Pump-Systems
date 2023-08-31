#coding: latin-1
#-------------------------------------------------------------------------------
# Exemple 6.1
#
from numpy import *
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse

fct = 'h_pompe_1510_10'
g = 9.81
#fluide = 'water'
fluide = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 40 % volume
patm = 101.325*1000.0
#
rho = 1024
nu = 3.28e-6
Cp = 3939
p1510 = array([-0.00403, -0.0057, 0.031,31.8])
rend_1510 = array([-0.00247, -0.818,13.6,14.33])

epsilon = 0
SDR = 11
D1 = 3
di1,do1 = sdr_pipe(D1,SDR)
A1 = pi*di1**2/4
L1 = 120.0
D2 = 1.25
di2,do2 = sdr_pipe(D2,SDR)
A2 = pi*di2**2/4
cas = 'c'
if cas == 'a':
    L2 = 435.0
    eff = 0.8
elif cas == 'b':
    L2 = 0.0
    eff = 0.8
elif cas == 'c':
    L2 = 0.0
    eff = 0.8*0.95
elif cas == 'd':
    L2 = 0.0
    eff = 0.8

nb = 10
V1 = A1*L1
m1 = V1*rho

def fct_h(q):
    if q == 0:
        return 0
    else:
        u1 = q/A1
        Re1 = di1*u1/nu
        ed1 = epsilon/di1
        f1 = Colebrook(Re1,ed1)
        h1 = f1*(L1/di1)/(A1**2*2*g)*q**2 # Dp en metres
        q2 = q/nb
        u2 = q2/A2
        Re2 = di2*u2/nu
        ed2 = epsilon/di2
        f2 = Colebrook(Re2,ed2)
        h1 = f1*(L1/di1)/(A1**2*2*g)*q**2 # Dp en metres
        h2 = f2*(L2/di2)/(A2**2*2*g)*q2**2 # Dp en metres
        hs = h1 + h2
        return hs
def fct2(q):
    hs = fct_h(q)
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp = h_pompe_1510(ql)
    y = hs - hp
    return y

def fct2(q):
    hs = fct_h(q)
    ql = q*1000.0   # q  m3/s,ql  litre/s
    hp = np.polyval(p1510,ql)
    y = hs - hp
    return y

p1510n = np.zeros(4)

def fct3(x,qx):
    hs = fct_h(qx)
    for i in range(0,4):
        p1510n[i] = p1510[i]*x**(i-1)
    ql = qx*1000.0   # q  m3/s,ql  litre/s
    hp = np.polyval(p1510n,ql)
    y = hs - hp
    return y

q = arange(0,160,2)  # gpm
ql = ls_gpm(q)
h1 = np.polyval(p1510,ql)
h1f = ft_m(h1)
if cas == 'a':
    qn = newton(fct2,15)
    mp = rho*qn
    qln = qn*1000       # debit en litre/s
    gpm = gpm_ls(qln)
    h_totn =  np.polyval(p1510,qln)
    rend = np.polyval(rend_1510,qln)
    rend = rend/100
    h_ft = ft_m(h_totn)
    print('flow (l/s) = ',qln,'gpm = ',gpm,' h = ',h_totn,' m')
    u1 = qn/A1
    Re1 = di1*u1/nu
    ed1 = epsilon/di1
    f1 = Colebrook(Re1,ed1)
    q2 = qn/nb
    u2 = q2/A2
    Re2 = di2*u2/nu
    ed2 = epsilon/di2
    f2 = Colebrook(Re2,ed2)
    hx1 = f1*(L1/di1)/(A1**2*2*g)*qn**2 # Dp en metres
    hx2 = f2*(L2/di2)/(A2**2*2*g)*q2**2 # Dp en metres
    hxs = hx1 + hx2
    W = mp*g*h_totn/rend/eff
    print('h (system)= ',hxs,' m')
    print('h (pump)= ',h_totn,' m')
    print('W  = ',W,' W')
elif cas == 'b':
    qn = newton(fct2,15)
    mp = rho*qn
    qln = qn*1000       # debit en litre/s
    gpm = gpm_ls(qln)
    h_totn =  np.polyval(p1510,qln)
    rend = np.polyval(rend_1510,qln)
    rend = rend/100
    h_ft = ft_m(h_totn)
    print('flow (l/s) = ',qln,'gpm = ',gpm,' h = ',h_totn,' m')
    u1 = qn/A1
    Re1 = di1*u1/nu
    ed1 = epsilon/di1
    f1 = Colebrook(Re1,ed1)
    hxs = f1*(L1/di1)/(A1**2*2*g)*qn**2 # Dp en metres
    print('h (system)= ',hxs,' m')
    print('h (pump)= ',h_totn,' m')
    W = mp*g*h_totn/rend/eff
    print('W  = ',W,' W')
elif cas == 'c':
    q_needed = 9.26/1000   # m3/s
    xx = newton(fct3,.5,args=(q_needed,))
    for i in range(0,4):
        p1510n[i] = p1510[i]*xx**(i-1)
    mp = rho*q_needed
    qln = q_needed*1000
    gpm = gpm_ls(qln)
    h_totn =  np.polyval(p1510n,qln)
    rend = np.polyval(rend_1510,qln/xx)
    rend = rend/100
    h_ft = ft_m(h_totn)
    print('flow (l/s) = ',qln,'gpm = ',gpm,' h = ',h_totn,' m')
    u1 =  q_needed/A1
    Re1 = di1*u1/nu
    ed1 = epsilon/di1
    f1 = Colebrook(Re1,ed1)
    hxs = f1*(L1/di1)/(A1**2*2*g)*q_needed**2 # Dp en metres
    print('h (system)= ',hxs,' m')
    print('h (pump)= ',h_totn,' m')
    W = mp*g*h_totn/rend/eff
    print('W  = ',W,' W')
elif cas == 'd':
    qn = newton(fct2,15)
    mp = rho*qn
    qln = qn*1000       # debit en litre/s
    gpm = gpm_ls(qln)
    h_totn =  np.polyval(p1510,qln)
    rend = np.polyval(rend_1510,qln)
    rend = rend/100
    h_ft = ft_m(h_totn)
    print('flow (l/s) = ',qln,'gpm = ',gpm,' h = ',h_totn,' m')
    u1 = qn/A1
    Re1 = di1*u1/nu
    ed1 = epsilon/di1
    f1 = Colebrook(Re1,ed1)
    hxs = f1*(L1/di1)/(A1**2*2*g)*qn**2 # Dp en metres
    print('h (system)= ',hxs,' m')
    print('h (pump)= ',h_totn,' m')
    W = mp*g*h_totn/rend/eff
    print('W  = ',W,' W')
    temps = 2*24*3600
    DT = W*temps/(m1*Cp)
    UA = 0.3*L1
    tau = m1*Cp/UA
    DT2 = W/UA*(1-exp(-temps/tau))
    print('DT = ',DT,DT2)

# plot system curve
#
n2 = len(ql)
h2 = zeros(n2)
for i in range(0,n2):
    h2[i] = fct_h(ql[i]/1000)
h2f = ft_m(h2)
plot(q,h1,q,h2)
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
rr = 4
rrn = 4
dx = rrn/2
dy = rrn*ratio/2
circle1 = Ellipse((gpm,h_ft), rr,rr*ratio, color='r')
ax.add_artist(circle1)
circle1.set_edgecolor('k')
xticks(fontsize=14)
yticks(fontsize=14)
show()
