#coding: latin-1
# Exemple 4.4
from geothermal_mod import *
from properties_mod import *
from numpy import *
from  matplotlib.pyplot import *
aljr = 0.10 # m2/jr
alhr = aljr/24.0
rb = 0.06
ks = 1.9
To = 10.0
#To = 0
q1 = -1000.0
q2 = 800.0
q3 = 500.0
q4 = -12000.0
Rb = 0.095
CCf = 2000
na = 20
Tfo = 25
Tfi = Tfo - q4/(CCf)
Tf = (Tfi+Tfo)/2.0
#Tf = -1.754
t1 = na*365
t2 = t1+60.0
t3 = t2+10.0
tf = t3+2.0
Fof = aljr*tf/rb**2
Fo1 = aljr*t1/rb**2
Fo2 = aljr*t2/rb**2
Fo3 = aljr*t3/rb**2
gam = 0.5772
G1a = g_function_lsi(Fof)
G2a = g_function_lsi(Fof-Fo1)
G3a = g_function_lsi(Fof-Fo2)
G4a = g_function_lsi(Fof-Fo3)
def gg(X):
    Y =(log(4*X)-0.5772)/(4*pi)
    return Y
def II(X):
    Y =(-log(X*X)-0.5772)/(2)
    return Y
G1 = (log(4*Fof)-gam)/(4*pi)
G2 = (log(4*(Fof-Fo1))-gam)/(4*pi)
G3 = (log(4*(Fof-Fo2))-gam)/(4*pi)
G4 = (log(4*(Fof-Fo3))-gam)/(4*pi)
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
q1p = q1/L
q2p = q2/L
q3p = q3/L
q4p = q4/L
Tb = To - (q1p*G1 + (q2p-q1p)*G2+ (q3p-q2p)*G3+ (q4p-q3p)*G4)/ks
Tfn = Tb - q4p*Rb
Tfon =  Tfn +  q4/(2*CCf)
print ('L = ',L)
print(q1*G1/ks/(To-Tf))
print((q2-q1)*G2/ks /(To-Tf))
print((q3-q2)*G3/ks/(To-Tf))
print((q4-q3)*G4/ks/(To-Tf))
print ('Tfo = ',Tfon)
print(Xf,X2,X3,X4)
print(G1,G2,G3,G4)
print(G1b,G2b,G3b,G4b)

d = 8
Fofb = aljr*tf/d**2
Fo1b = aljr*t1/d**2
Fo2b = aljr*t2/d**2
Fo3b = aljr*t3/d**2
G1bb = gg(Fofb)
G2bb = g_function_lsi(Fofb-Fo1b)
G3bb = gg(Fofb-Fo2b)
G4bb = gg(Fofb-Fo3b)
DT2a = -(q1p*(G1bb-G2bb))/ks
print (DT2a)
DT2b = -q1p*log(tf/(tf-t1))/(4*pi*ks)
print (DT2b)
X1b  = d/(2.0*sqrt(aljr*tf))
X2b  = d/(2.0*sqrt(aljr*(tf-t1)))
X3b  = d/(2.0*sqrt(aljr*(tf-t2)))
X4b  = d/(2.0*sqrt(aljr*(tf-t3)))
I1  = II(X1b)
I2  = fct_I(X2b)
I3  = fct_I(X3b)
I4  = fct_I(X4b)
DT1 = -q1p/(2*pi*ks)*(I1)
DT2 = -(q2p-q1p)/(2*pi*ks)*I2
DT3 = -(q3p-q2p)/(2*pi*ks)*I3
DT4 = -(q4p-q3p)/(2*pi*ks)*I3
DT = DT1+DT2+DT3+DT4       # Caclul de la perturbaion totale
print  (X1b,X2b,X3b,X4b)
print  (DT1,DT2,DT3,DT4,DT)

Ra = 0.1
eta = L/(CCf*sqrt(Rb*Ra))
Rbs = Rb*eta/tanh(eta)
Tfn = Tb - q4p*Rbs
Tfon = Tfn + q4/(2*CCf)
print ('eta = ',eta)
print ('Rbs = ',Rbs)
print ('Tfon = ',Tfon)
