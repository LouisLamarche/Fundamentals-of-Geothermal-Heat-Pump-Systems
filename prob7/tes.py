#
# exemple de la methode de la ligne source pour évaluer la conductivité du sol
#
from numpy import *
from geothermal_md import *
from conversion_md import *
from hydraulic_md import *
from scipy.optimize import newton
Patm = 101.325*1000.0
#
#
#
# fluide
fluide = 'Water'
Trefk = 290
mu = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
rho = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Pr = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
nu = mu/rho
qv = array([3.9,7,13,19,35,100])
qm = array([3.7,5.7,9,12,18,40])
qmax = ls_gpm(qv)
qmin = ls_gpm(qm)
tu =array([0.75,1,1.25,1.5,2,3])
n = len(qv)
uv = zeros(n)
gmax = zeros(n)
gmin = zeros(n)
SDR = 11
L = m_ft(100)
g = 9.81
epp = 0.00/1000
def calcul_gpm(x):
    ux = x/A1
    Rex = ux*d1/nu
    f = Colebrook(Rex,ed)
    h = f*(L/d1)*ux**2/(2*g)
    hp = ft_m(h)
    y = hp-4
    return(y)

for i in range(0,n):
    d1,d2 = sdr_pipe(tu[i],SDR)     # choix des tuyaux SDR-11 1.25 po nominal
    A1 = pi*d1**2/4.0
    Q1 = mcs_gpm(qv[i])
    u1 = Q1/A1
    Re = u1*d1/nu
    ed = epp/d1
    deb = newton(calcul_gpm,Q1)
    gmax[i] = gpm_mcs(deb)
    u2 = m_ft(2)
    Q2 = u2*A1
    gmin[i] = gpm_mcs(Q2)
