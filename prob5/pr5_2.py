#coding: utf-8
from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
#
data = loadtxt("pr2a.txt")  # charges en kWatts
ta = data[:,0]    # températures extérieures
Tout_tra = data[:,1]    # températures extérieures
q_tra= data[:,2]    # températures extérieures
data = loadtxt("pr2b.txt")  # charges en kWatts
tb = data[:,0]    # températures extérieures
Tout_trb = data[:,1]  + 10   # températures extérieures
q_trb= data[:,2]    # températures extérieures
data = loadtxt("pr2ai.txt")  # charges en kWatts
tai = data[:,0]    # températures extérieures  avec                             short-circuit
Tout_trai = data[:,1]    # températures extérieures
q_trai= data[:,2]    # températures extérieures
data = loadtxt("pr2bi.txt")  # charges en kWatts
tbi = data[:,0]    # températures extérieures
Tout_trbi = data[:,1]  + 10  # températures extérieures
q_trbi = data[:,2]    # températures extérieures

# ground
als = 0.1 # m2/jr
alhr = als/24.0
rb = 0.055
ks = 2.5
To = 10.0
H = 200.0
# puits
d1,d2 = sdr_pipe(1,11)
ro= d2/2.0
ri = d1/2.0
xc = 0.02
kg = 1.2
kp = 0.4
mpa = 2000.0/3600.0
mpb = 1000.0/3600.0
Cp = 4190.0
rho = 1000.0
q1p = -35.0
q2p = 20.0
q3p = -60.0
# time intermédiaires
t1 = 120.0
t2 = 350.0
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
tf1 = 250
tf2 = 400

cas = 'a'
if cas == 'a':
    mp =  mpa
    t = ta
    Tout_tr =Tout_tra
    Tout_tri =Tout_trai
else:
    t = tb
    mp =  mpb
    Tout_tr =Tout_trb
    Tout_tri =Tout_trbi
#
#
#
def calcul_Tb(tf):
    Fo = alhr*tf/rb**2
    if tf < t1:
        DT1 = -q1p*G_function(Fo)/ks
        qp = q1p
    elif tf < t2:
        DT1 = -(q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1))/ks
        qp = q2p
    else:
        DT1 = -(q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1) + (q3p - q2p)*G_function(Fo-Fo2))/ks
        qp = q3p
    return DT1,qp
#
# Calcul de Rb
#
nu = 8.54e-7
Ac = pi*d1**2/4.0
um = mp/(rho*Ac)
Re = um*d1/nu
Pr = 5.86
kf = 0.61
if (Re>2300.0):
        # Gnielienski
        f = (0.79*log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    print(' laminar')
hf = (Nud*kf)/(d1)
rconv = 1/(pi*d1*hf)
rcond = log(ro/ri)/(2*pi*kp)
Rp = rcond + rconv
print('R conduction = ',rcond)
print('R convection = ',rconv)
print('R conduite total = ',Rp)
Rbs = Rb_Sharqawi(kg,rb,ro,xc)
Rbh,Rah = Rb_linesource(kg,ks,rb,ro,xc)
Rbs = Rbs + Rp/2.0
Rbh = Rbh + Rp/2.0
Rah = Rah + 2.0*Rp
print('Rb Hellstrom (xc = 3.5) =  ',Rbh)
print('Rb Sharq = (xc = 3.5) =  ',Rbs)
#
# Calcul de Tb t = 250
#
#b
#
gamh = H/(mp*Cp*Rbh)
xh = gamh/2
Rbhb = Rbh*xh/tanh(xh)
gams = H/(mp*Cp*Rbs)
xs = gams/2
Rbsb = Rbs*xs/tanh(xs)
#c
#
eta = H/(mp*Cp*sqrt(Rah*Rbh))
Rbhc = Rbh*eta/tanh(eta)
Rbhca = Rbh*(1 + eta**2/3)
print(Rbh)
print(Rbsb)
print(Rbhc)
R12 = 4*Rah*Rbh/(4*Rbh-Rah)
print(R12)
fl_pl = 1
if fl_pl:
    nn = len(t)
    Tbv = zeros(nn)
    qpv = zeros(nn)
    for i in range(0,nn):
        DT,qp = calcul_Tb(t[i])
        Tbv[i] = To +  DT
        qpv[i] = qp
    Tfv = Tbv - qpv*Rbhb
    Tfov = Tfv + qpv*H/(2.0*mp*Cp)
    Tfvi = Tbv - qpv*Rbhc
    Tfovi = Tfvi + qpv*H/(2.0*mp*Cp)
    figure(1)
    plot(t,Tfov,t,Tout_tr)
    legend(('Tout calc','Tout TRSYS '))
    figure(2)
    plot(t,Tfovi,t,Tout_tri)
    legend(('Tin calc',' Tfi TRSYS'))
    show()
DT1a,qpa = calcul_Tb(tf1)
Tb1a = To + DT1a  # Calcul de Tb avec la  source cylindrique
print('')
print('Tb  (t = 250 heures) =  ',Tb1a)
#
Tfia =  Tb1a - qpa*Rbh
Tfiia =  Tb1a - qpa*Rbs
Tfoi = Tfia + qpa*H/(2.0*mp*Cp)
Tfoii = Tfiia + qpa*H/(2.0*mp*Cp)
Tfib =  Tb1a - qpa*Rbhb
Tfiib =  Tb1a - qpa*Rbsb
Tfoiii = Tfib + qpa*H/(2.0*mp*Cp)
Tfoiv = Tfiib + qpa*H/(2.0*mp*Cp)
Tfic =  Tb1a - qpa*Rbhc
Tfov = Tfic + qpa*H/(2.0*mp*Cp)
print ('\t Tout(250 heures) a i) = ' , Tfoi)
print ('\t Tout(250 heures) a ii) = ' , Tfoii)
print ('\t Tout(250 heures) a iii) = ' , Tfoiii)
print ('\t Tout(250 heures) a iv) = ' , Tfoiv)
print ('\t Tout(250 heures) a v) = ' , Tfov)
print ('\t Tout(250 heures) Trnsys = ' ,Tout_tr[250])
print ('\t Tout(250 heures) Trnsys inter = ' ,Tout_tri[250])

DT1b,qpb = calcul_Tb(tf2)
Tb1b = To + DT1b  # Calcul de Tb avec la  source cylindrique
print('')
print('Tb (t = 400 heures) =  ',Tb1b)
Tfia =  Tb1b - qpb*Rbh
Tfiia =  Tb1b - qpb*Rbs
Tfoi = Tfia + qpb*H/(2.0*mp*Cp)
Tfoii = Tfiia + qpb*H/(2.0*mp*Cp)
Tfib =  Tb1b - qpb*Rbhb
Tfiib =  Tb1b - qpb*Rbsb
Tfoiii = Tfib + qpb*H/(2.0*mp*Cp)
Tfoiv = Tfiib + qpb*H/(2.0*mp*Cp)
Tfic =  Tb1b - qpb*Rbhc
Tfov = Tfic + qpb*H/(2.0*mp*Cp)
print ('\t Tout(400 heures) a i) = ' , Tfoi)
print ('\t Tout(400 heures) a ii) = ' , Tfoii)
print ('\t Tout(400 heures) a iii) = ' , Tfoiii)
print ('\t Tout(400 heures) a iv) = ' , Tfoiv)
print ('\t Tout(400 heures) a v) = ' , Tfov)
print ('\t Tout(400 heures) Trnsys = ' ,Tout_tr[400])
print ('\t Tout(400 heures) Trnsys = ' ,Tout_tri[400])
