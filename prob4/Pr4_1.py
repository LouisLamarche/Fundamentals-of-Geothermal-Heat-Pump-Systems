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
tf = 150 -ti
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
Fo = alhr*tf/rb**2
def G_approx(Fo):
    return 1/(4*pi)*(np.log(4*Fo) - gam)
def I_approx(X):
    return -1/2*(2*np.log(X) + gam)
#
#a
#
DTi = -(q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1)+(q3p-q2p)*G_function_ils(Fo-Fo2))/ks
DT = -(q1p*G_approx(Fo)+(q2p-q1p)*G_approx(Fo-Fo1)+(q3p-q2p)*G_approx(Fo-Fo2))/ks
qp = q3p
#
# or
#
Xf = rb/(2*sqrt(alhr*tf))
X1 = rb/(2*sqrt(alhr*(tf-t1)))
X2 = rb/(2*sqrt(alhr*(tf-t2)))
DTiX = -(q1p*I_function(Xf)+(q2p-q1p)*I_function(X1)+(q3p-q2p)*I_function(X2))/(2*pi*ks)
DTX = -(q1p*I_approx(Xf)+(q2p-q1p)*I_approx(X1)+(q3p-q2p)*I_approx(X2))/(2*pi*ks)
print('DT a) LSI (t = 150 heures) =  ',DT,DTX)

#
# b
#
qagg = (q1p*t1 + q2p*(t2 - t1))/t2
DTib = -(qagg*G_function_ils(Fo)+(q3p-qagg)*G_function_ils(Fo-Fo2))/ks
DTb = -(qagg*G_approx(Fo)+(q3p-qagg)*G_approx(Fo-Fo2))/ks
# or
DTiXb = -(qagg*I_function(Xf)+(q3p-qagg)*I_function(X2))/(2*pi*ks)
DTXb = -(qagg*I_approx(Xf)+(q3p-qagg)*I_approx(X2))/(2*pi*ks)
print('DT b) LSI (t = 150 heures) =  ',DTb,DTXb)
