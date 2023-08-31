#coding: latin-1
# Exemple 4.4
from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
aljr = 0.10 # m2/jr
alhr = aljr/24.0
rb = 0.06
ks = 1.9
To = 10.0
L = 200
q1 = -1000.0
q2 = 800.0
q3 = 500.0
q4 = -12000.0
na = 20
t1 = na*365
t2 = t1+60.0
t3 = t2+10.0
tf = t3+2.0
Fof = aljr*tf/rb**2
Fo1 = aljr*t1/rb**2
Fo2 = aljr*t2/rb**2
Fo3 = aljr*t3/rb**2
gam = euler_gamma
G1a = G_function_ils(Fof)
G2a = G_function_ils(Fof-Fo1)
G3a = G_function_ils(Fof-Fo2)
G4a = G_function_ils(Fof-Fo3)
def gg(X):
    Y =(log(4*X)-gam)/(4*pi)
    return Y
def II(X):
    Y =(-log(X*X)-gam)/(2)
    return Y
G1b = gg(Fof)
G2b = gg(Fof-Fo1)
G3b = gg(Fof-Fo2)
G4b = gg(Fof-Fo3)
Xf = rb/(2*sqrt(aljr*tf))
X2 = rb/(2*sqrt(aljr*(tf-t1)))
X3 = rb/(2*sqrt(aljr*(tf-t2)))
X4 = rb/(2*sqrt(aljr*(tf-t3)))
G1c = II(Xf)/(2*pi)
G2c = II(X2)/(2*pi)
G3c = II(X3)/(2*pi)
G4c = II(X4)/(2*pi)
G1d = G_function(Fof)
G2d = G_function(Fof-Fo1)
G3d = G_function(Fof-Fo2)
G4d = G_function(Fof-Fo3)
q1p = q1/L
q2p = q2/L
q3p = q3/L
q4p = q4/L
Tba = To - (q1p*G1a + (q2p-q1p)*G2a+ (q3p-q2p)*G3a+ (q4p-q3p)*G4a)/ks
Tbb = To - (q1p*G1b + (q2p-q1p)*G2b+ (q3p-q2p)*G3b+ (q4p-q3p)*G4b)/ks
Tbc = To - (q1p*G1c + (q2p-q1p)*G2c+ (q3p-q2p)*G3c+ (q4p-q3p)*G4c)/ks
Tbd = To - (q1p*G1d + (q2p-q1p)*G2d+ (q3p-q2p)*G3d+ (q4p-q3p)*G4d)/ks
print ('Tb = ',Tba,Tbb,Tbc,Tbd)
qagga = (q2p*(t2-t1) + q3p*(t3 - t2))/(t3-t1)
Tbdna = To - (q1p*G1d + (qagga-q1p)*G2d + (q4p-qagga)*G4d)/ks
qagg = (q1p*t1 + q2p*(t2 - t1) + q3p*(t3 - t2))/t3
Tbdn = To - (qagg*G1d + (q4p-qagg)*G4d)/ks
print ('Tb = ',Tbd,Tbdna,Tbdn)


