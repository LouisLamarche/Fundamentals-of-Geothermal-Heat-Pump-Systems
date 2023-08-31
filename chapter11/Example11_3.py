#coding: utf-8
#
# exemple 11-3
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
#
M = np.loadtxt("..\\data\\pump.txt")
tm = M[:,0]     # temps en minutes
t = tm      # t en min
st = M[:,1]   # rabattement en metres
s1 = M[:,2]   # rabattement en metres
s2 = M[:,3]   # rabattement en metres
q  = 0.008     #  m3/min
q = q*60
rw = 0.30/2.0
r1 = 5.0
r2 = 25
nt = len(st)
n1 = 18
#
# On boucle en prenant les N-i derniers points pour faire la régression
#
xtm = np.log(tm)
xt = np.log(t)
x = np.log(t[n1:nt-1])
y = st[n1:nt-1]
y1 = s1[n1:nt-1]
y2 = s2[n1:nt-1]
p = np.polyfit(x,y,1)
p1 = np.polyfit(x,y1,1)
p2 = np.polyfit(x,y2,1)

m0 = p[0]
v0 = p[1]
m1 = p1[0]
v1 = p1[1]
m2 = p2[0]
v2 = p2[1]
to = np.exp(-v0/m0)
to1 = np.exp(-v1/m1)
to2 = np.exp(-v2/m2)
Ts = q/(4.0*pi*m0)
Ts1 = q/(4.0*pi*m1)
Ts2 = q/(4.0*pi*m2)
s2s = s2[nt-1]
s1s = s1[nt-1]
Tss = q/(2*pi)/(s2s-s1s)*np.log(r1/r2)
print ('T (regime permanent ) = ',Tss,' m2/min')
gam = 0.5772157
sws1 = s1s + q/(2*pi*Tss)*np.log(r1/rw)
sws2 = s2s + q/(2*pi*Tss)*np.log(r2/rw)
print ('sw (regime permanent ) = ',sws1,' m')
ro = np.exp(sws1*2*pi*Tss/q)*rw
print ('ro (regime permanent ) = ',ro,' m')
ro = np.exp(s1s*2*pi*Tss/q)*r1
print ('ro (regime permanent ) = ',ro,' m')
ro = np.exp(sws2*2*pi*Tss/q)*rw
beta = np.log(ro/rw)/(2*pi*Tss)
print ('beta = ',beta,'min/m2')
sts = st[nt-1]
Css = (sts - sws1)/q**2
print ('C (régime permanent) = ',Css,' min2/m5')
Css2 = (sts - beta*q)/q**2
So = 4*Ts*to*np.exp(-gam)/rw**2     # aucun sens
S1 = 4*Ts1*to1*np.exp(-gam)/r1**2
S2 = 4*Ts2*to2*np.exp(-gam)/r2**2
C = (gam-np.log(4*to*Ts/(S1*rw**2)))/(4*pi*Ts*q)
C2 = (gam-np.log(4*to*Ts/(S2*rw**2)))/(4*pi*Ts*q)
print ('C (régime transitoire) = ',C,' min2/m5')
print ('C (régime transitoire) = ',C2,' min2/m5')
st2 = np.polyval(p,xt)
s12 = np.polyval(p1,xt)
s22 = np.polyval(p2,xt)
print ('T (puits o) = ',Ts,' m2/min')
print ('T (puits 1) = ',Ts1,' m2/min')
print ('T (puits 2) = ',Ts2,' m2/min')

print ('S (puits o) = ',So)
print ('S (puits 1) = ',S1)
print ('S (puits 2) = ',S2)

print ('al (puits o) = ',Ts/So)
print ('al (puits 1) = ',Ts1/S1)
print ('al (puits 2) = ',Ts2/S2)

alo = Ts/So
al1 = Ts1/S1
al2 = Ts2/S2
t11 = t[n1]
Foo = al1*t11/rw**2 ; print ( 'Foo = ',Foo)
Fo1 = al1*t11/r1**2 ; print ( 'Foo = ',Fo1)
Fo2 = al1*t11/r2**2 ; print ( 'Foo = ',Fo2)
uoo = 1/(4*Foo) ; print ( 'uoo = ',uoo)
uo1 = 1/(4*Fo1) ; print ( 'uoo = ',uo1)
uo2 = 1/(4*Fo2); print ( 'uoo = ',uo2)
plot(xtm,st,'*',xtm,st2,xtm,s1,'*',xtm,s12,xtm,s2,'*',xtm,s22)
grid()
ylabel('s(m)',fontsize = 15,fontweight ='bold')
xlabel('log(t)',fontsize = 15,fontweight ='bold')
show()



