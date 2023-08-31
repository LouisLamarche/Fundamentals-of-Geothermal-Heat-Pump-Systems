#coding: latin-1
# Exemple 4.4
from geothermal_md import *
import numpy as np
from matplotlib.pyplot import *
#
# data
#
aljr = 0.10 # m2/jr
alhr = aljr/24.0
rb = 0.06
dtt = alhr*1/rb**2
ks = 2.5
To = 10.0
d = 6
nx = 3
ny = 2
nt = nx*ny
zt = np.zeros([nt,2])
k = 0
H = 100
Bt = d/H
rr = rb/H
dth = 8*dtt*rr**2
for i1  in range(0,nx):
    x = i1*Bt
    for i2 in range(0,ny):
        y = i2*Bt
        zt[k] = [x,y]
        k = k+1
#
#############################################################################################
#############################################################################################
#
# calcul de la fonction g analytique associée au champ de capteur choisi
#
#
#g[i] = compute_g_function(zt,tt[i],rbb = rr)    # SI on n'ecrit pas D/H , la fonction prend 0.04 par défaut ( Comme Eskilson)

q1h = 15000.0
q2h = -14000.0
q3h = 13000.0
cop = 3
q1 = q1h*(cop-1)/cop
q2 = q2h*(cop+1)/cop
q3 = q3h*(cop-1)/cop
Rb = 0.08
CCf = 2850
Tf = -2

#
# modele sans aggregation
#
flag = 0
nt = 48
if flag == 0:
    q1r = np.random.randint(2, size=nt)
    q2r = np.random.randint(2, size=nt)
    q3r = np.random.randint(2, size=nt)
    q = q1*q1r + q2*q2r + q3*q3r
    np.savetxt('load.txt',q)
else:
    q = np.loadtxt('load.txt')
n2  = int(nt/4)
q1n = np.zeros(nt)
q2n = np.zeros(nt)
q3n = np.zeros(nt)
for i in range(0,n2):
    i1 = i*4
    i2 = i1+4
    q1m = q1*np.mean(q1r[i1:i2])
    if q1m> q1/4:
        q1n[i1:i2] = q1
    q2m = q2*np.mean(q2r[i1:i2])
    if q2m < q2/4:
        q2n[i1:i2] = q2
    q3m = q3*np.mean(q3r[i1:i2])
    if q3m> q3/4:
        q3n[i1:i2] = q3
t = np.arange(0,nt)
plot(t,q3*q3r,t,q3n)
show()
exit()
#q[0:nt-1:2] = q1
#q[nt-1] = q1
fctG = G_function_ils
K = np.zeros(nt)
L =  0
hr = 0
#Fom = dtt*(hr+1)
Fom = dth*(hr+1)
den = 2*pi
#gM= fctG(Fom)
gM = compute_g_function(zt,Fom,rbb = rr)    # SI on n'ecrit pas D/H , la fonction prend 0.04 par défaut ( Comme Eskilson)
K[hr] = gM
for hr in range(1,nt):
#    Fop = dtt*(hr+1)
    Fop = dth*(hr+1)
#    gK= fctG(Fop)
    gK = compute_g_function(zt,Fop,rbb = rr)    # SI on n'ecrit pas D/H , la fonction prend 0.04 par défaut ( Comme Eskilson)
    K[hr] = gK-gM
    gM = gK
for hr in range(0,nt):
    s = 0
    for j in range(0,hr+1):
        print(hr-j,K[hr-j])
        s=s+q[j]*K[hr-j]
    SCT = s/(den*ks)
    Ln = (SCT + q[hr]*Rb)/(To-Tf)
    L = max(L,Ln)


print ('L = ',L)
#exit()
Tw1 = np.zeros(nt)
for hr in range(0,nt):
    s = 0
    for j in range(0,hr+1):
        s=s+q[j]*K[hr-j]
    Tw1[hr] = To - s/(L*den*ks) - q[hr]*Rb/L
plot(Tw1)
show()
exit()

##print ('fin du calcul  ',toc)
##print(' fin du modele sans aggrégation 2')

def gg(X):
    Y =(np.log(4*X)-0.5772)/(4*pi)
    return Y
def II(X):
    Y =(-np.log(X*X)-0.5772)/(2)
    return Y
G1a = G_function(Fof)
G2a = G_function(Fof-Fo1)
G3a = G_function(Fof-Fo2)
G4a = G_function(Fof-Fo3)
G1b = gg(Fof)
G2b = gg(Fof-Fo1)
G3b = gg(Fof-Fo2)
G4b = gg(Fof-Fo3)
Xf = rb/(2*np.sqrt(aljr*tf))
X2 = rb/(2*np.sqrt(aljr*(tf-t1)))
X3 = rb/(2*np.sqrt(aljr*(tf-t2)))
X4 = rb/(2*np.sqrt(aljr*(tf-t3)))
G1 = II(Xf)/(2*pi)
G2 = II(X2)/(2*pi)
G3 = II(X3)/(2*pi)
G4 = II(X4)/(2*pi)
SCT = (q1*G1 + (q2-q1)*G2+ (q3-q2)*G3+ (q4-q3)*G4)/ks
L = (SCT + q4*Rb)/(To-Tf)
print ('L = ',L)
#print(Xf,X2,X3,X4)
#print(G1,G2,G3,G4)
#print(G1a,G2a,G3a,G4a)
#print(G1b,G2b,G3b,G4b)
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
X1b  = d/(2.0*np.sqrt(aljr*tf))
I1  = II(X1b)
DT1 = -qp1/(2*pi*ks)*(I1)
Fob = aljr*tf/d**2
G1 = gg(Fob)
DT2 = -qp1/ks*G1
#print  (X1b,Fob)
print  ('b) Delt T = ',DT1,DT2)
#
# c)
#
Ra = 0.1
eta = L/(CCf*np.sqrt(Rb*Ra))
Rbs = Rb*eta/np.tanh(eta)
Tfn = Tb - qp4*Rbs
Tfonn = Tfn + q4/(2*CCf)
print ('eta = ',eta)
print ('Rbs = ',Rbs)
print ('Tfon = ',Tfonn)
