#coding: utf-8
#
# Example 4.4 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
Rp = 0.08
ro = 0.015
Rpp = 0.01
h = 1/Rpp
Rp = 1/(2*pi*ro*h)
ks = 2.5   # soil
kg = ks   # grout
rhos = 1000
rb = 0.08
xc = 0.042
Cp_s = 2500
als = ks/(rhos*Cp_s)
alhr = als*3600
alj = alhr*24
H = 150
q = -7500
qp = q/H
To = -2.5
To = 0
data = np.loadtxt('Tp_hart_couvillon.txt')
t1 = data[:,0]
Tp = data[:,1]
data = np.loadtxt('Tf_hart_couvillon.txt')
t2 = data[:,0]
Tf = data[:,1]
data = np.loadtxt('Tb_hart_couvillon.txt')
t3 = data[:,0]
Tb = data[:,1]
#
# line source
#
nn= len(t1)
dt = 10
t = t1[1:nn+dt:dt]
Tbn = Tb[1:nn+dt:dt]
Tpn = Tp[1:nn+dt:dt]
Tfn = Tf[1:nn+dt:dt]
nb = len(t)
Tbnn = np.zeros(nb)
Tpnn = np.zeros(nb)
Tfnn = np.zeros(nb)
Tpnh1 = np.zeros(nb)
Tpnh2 = np.zeros(nb)
z = np.array([xc,-xc])  #
J = 10
Rg,Rga = Rb_multipole(kg,ks,rb,ro,0,J,z)
Rg,Rga = Rb_linesource(kg,ks,rb,ro,xc)
for i in range(nb):
    tt = t[i]
    Fo =  alhr*tt/rb**2
    Tb = To - qp/ks*G_function(Fo)
    Tbnn[i] = Tb
    sigma = (kg-ks)/(kg+ks)
    Rb1 = Rg + Rp/2.0
    Ra1 = Rga + 2*Rp
    Tp = Tb - qp*Rg
    Tpnn[i] = Tp
    Tf = Tb - qp*Rb1
    Tfnn[i] = Tf
    Tf2 = Tp - qp*Rp/2
    Foo = alhr*tt/ro**2
    G1,r =  g_function_hart_couv(Foo)
    sd = 2*xc
    Fo2 = alhr*tt/sd**2
    Fo3 = alhr*tt/xc**2
    G2,r =  g_function_hart_couv(Fo2)
    G3,r =  g_function_hart_couv(Fo3)
    LMLS1 = 1 + G2/G1
    DTm = 0.5*G3/G1
    rinf = 4*np.sqrt(alhr*tt)
    y = 1/(4*np.sqrt(Fo2))
    the = np.arccos(y/2)
    F = DTm*(the -np.cos(the)*np.sin(the))/pi
    LMLS2 = 1/(1-F)
    LMLS = LMLS1*LMLS2
    qp2 = qp/2
    Tph1 = To - qp2/ks*G_function(Foo)*LMLS1
    Tph2 = To - qp2/ks*G_function(Foo)*LMLS
    Tpnh1[i] = Tph1
    Tpnh2[i] = Tph2
    Tfh1 = Tph1 - qp2*Rp
    Tfh2 = Tph2 - qp2*Rp
p1 = semilogx(t,Tpn,color = 'k',linestyle = '-',marker = 'x',markersize=11,label = 'Comsol ')
p2 = semilogx(t,Tpnn,color = 'k',linestyle = '-',marker = '+',markersize=11,label = 'Borehole resistance')
p3 = semilogx(t,Tpnh1,color = 'k',linestyle = '-',marker = 'o',markersize=11,label = 'Hart-Couvillon ($LMLS_{HE}$)')
p4 = semilogx(t,Tpnh2,color = 'k',linestyle = '-',marker = '1',markersize=11,label = 'Hart-Couvillon ($LMLS$)')
legend(fontsize = 14)
ax = gca()
grid(True,which='both')
#title('Coubes de pertes de charge  ')
fts = 16
ftst = 14
xlabel(' Fo',fontname='Times new Roman',fontsize = fts)
ylabel('G',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()


