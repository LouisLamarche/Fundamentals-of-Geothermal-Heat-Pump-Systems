#coding: latin-1
# Exemple 4.4
from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
#
# data
#
aljr = 0.10 # m2/jr
alhr = aljr/24.0
rb = 0.06
ks = 2.0
To = 7.5
d = 6
q1 = 1500.0
q2 = -2800.0
q3 = 500.0
q4 = 5000.0
Rb = 0.095
CCf = 2850
na = 20
dt1 = 60
dt2 = 30
dt3 = 1
Tfo = 0
Tfi = Tfo - q4/(CCf)
Tf = (Tfi+Tfo)/2.0
t1 = na*365
t2 = t1+dt1
t3 = t2+dt2
tf = t3+dt3
Fof = aljr*tf/rb**2
Fo1 = aljr*t1/rb**2
Fo2 = aljr*t2/rb**2
Fo3 = aljr*t3/rb**2
gam = 0.5772
def gg(X):
    Y =(log(4*X)-0.5772)/(4*pi)
    return Y
def II(X):
    Y =(-log(X*X)-0.5772)/(2)
    return Y
G1a = G_function(Fof)
G2a = G_function(Fof-Fo1)
G3a = G_function(Fof-Fo2)
G4a = G_function(Fof-Fo3)
G1b = gg(Fof)
G2b = gg(Fof-Fo1)
G3b = gg(Fof-Fo2)
G4b = gg(Fof-Fo3)
Xf = rb/(2*sqrt(aljr*tf))
X2 = rb/(2*sqrt(aljr*(tf-t1)))
X3 = rb/(2*sqrt(aljr*(tf-t2)))
X4 = rb/(2*sqrt(aljr*(tf-t3)))
G1 = II(Xf)/(2*pi)
G2 = II(X2)/(2*pi)
G3 = II(X3)/(2*pi)
G4 = II(X4)/(2*pi)
SCT = (q1*G1 + (q2-q1)*G2+ (q3-q2)*G3+ (q4-q3)*G4)/ks
L = (SCT + q4*Rb)/(To-Tf)
print ('L = ',L)
print(Xf,X2,X3,X4)
print(G1,G2,G3,G4)
print(G1a,G2a,G3a,G4a)
print(G1b,G2b,G3b,G4b)
#
# verification
#
qp1 = q1/L
qp2 = q2/L
qp3 = q3/L
qp4 = q4/L
H = L/2.0
Tb = To -(qp1*G1 + (qp2-qp1)*G2+ (qp3-qp2)*G3+ (qp4-qp3)*G4)/ks
Tf = Tb - qp4*Rb
Tfon = Tf + q4/(2*CCf)
print ('Tfo = ',Tfon)
#
# b)
#
qpm = (q1*t1 + q2*dt1 + q3*dt2 + q4*dt3)/tf
X1b  = d/(2.0*sqrt(aljr*tf))
I1  = II(X1b)
DT1 = -qp1/(2*pi*ks)*(I1)
Fob = aljr*tf/d**2
G1 = gg(Fob)
DT2 = -qp1/ks*G1
print  (X1b,Fob)
print  ('Delt T = ',DT1,DT2)
#
# c)
#
Ra = 0.1
eta = L/(CCf*sqrt(Rb*Ra))
Rbs = Rb*eta/tanh(eta)
Tfn = Tb - qp4*Rbs
Tfonn = Tfn + q4/(2*CCf)
print ('eta = ',eta)
print ('Rbs = ',Rbs)
print ('Tfon = ',Tfonn)
