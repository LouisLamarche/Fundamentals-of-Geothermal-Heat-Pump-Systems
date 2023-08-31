#coding: latin-1
from geothermal_md import *
from numpy import *
#

# sol
gam = np.euler_gamma   # euler constant
pi = np.pi
als = 0.088
alhr = als/24.0
ks = 2.25
q1p = 10
q2p = 40.0
q3p = 100.0
# valeurs des temps intermédiaires
ti = 7
t0 = 7 - ti
t1 = 50 - ti
t2 = 115 - ti
rb = 0.145
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
X1 = rb/(2*sqrt(alhr*t1))
X2 = rb/(2*sqrt(alhr*t2))
tf1 = 150 -ti
def G_approx(Fo):
    return 1/(4*pi)*(np.log(4*Fo) - gam)
def I_approx(X):
    return -1/2*(2*np.log(X) + gam)
#
#
def calcul_Tb(tf):
    Fo = alhr*tf/rb**2
    if tf < t1:
        DTi = -q1p*G_function_ils(Fo)/ks
        DT = -q1p*G_approx(Fo)/ks
        qp = q1p
    elif tf < t2:
        DTi = -(q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1))/ks
        DT = -(q1p*G_approx(Fo)+(q2p-q1p)*G_approx(Fo-Fo1))/ks
        qp = q2p
    else:
        DTi = -(q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1)+(q3p-q2p)*G_function_ils(Fo-Fo2))/ks
        DT = -(q1p*G_approx(Fo)+(q2p-q1p)*G_approx(Fo-Fo1)+(q3p-q2p)*G_approx(Fo-Fo2))/ks
        qp = q3p
    return DT,qp
def calcul_TbX(tf):
    Xf = rb/(2*sqrt(alhr*tf))
    X1 = rb/(2*sqrt(alhr*(tf-t1)))
    X2 = rb/(2*sqrt(alhr*(tf-t2)))
    if tf < t1:
        DTi = -q1p*I_function(Xf)/(2*pi*ks)
        DT = -q1p*I_approx(Fo)/(2*pi*ks)
        qp = q1p
    elif tf < t2:
        DTi = -(q1p*I_function(Xf)+(q2p-q1p)*I_function(X1))/(2*pi*ks)
        DT = -(q1p*I_approx(Xf)+(q2p-q1p)*I_approx(X1))/(2*pi*ks)
        qp = q2p
    else:
        DTi = -(q1p*I_function(Xf)+(q2p-q1p)*I_function(X1)+(q3p-q2p)*I_function(X2))/(2*pi*ks)
        DT = -(q1p*I_approx(Xf)+(q2p-q1p)*I_approx(X1)+(q3p-q2p)*I_approx(X2))/(2*pi*ks)
        qp = q3p
    return DT,qp
DT1,qpa = calcul_Tb(tf1)
#
# or
#
DT2,qpa = calcul_TbX(tf1)
print('DTb LSI (t = 150 heures) =  ',DT1,DT2)

