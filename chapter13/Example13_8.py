#coding: utf-8
#
# Exemple 13.8
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from time import *
# fichier de fonction g (Eskilson) tabuléees pour champ 2 x 2 pour b = 0.05,0.1,0.2,0.4,0.8
zo = 0
ro = 0.25/2
b = 0.5
H = 5
Ht = H/ro
bt = b/ro
Fo = 5
rp = 0.01
rpt = rp/ro
zot = zo/ro
ome = 2*pi/b
omet = ome*ro
z = H/2
phi = ome*z
phin = ome*H/4
Fov = np.linspace(5,200,20)
n = len(Fov)
g1 = np.zeros(n)
g3 = np.zeros(n)
g2 = np.zeros(n)
phi1 = 0
phi2 = ome*H
def g_integrand(x,Fo,rpt,Ht,zot,bt):
    rr = 1 + rpt*cos(x)
    zz = phi/omet + rpt*sin(x)
    y = G_spiral_pipe(Fo,rpt,x,Ht,zot,bt)
    return y
for i in range(0,n):
    Fo = Fov[i]
    g1[i] = G_spiral_pipe(Fo,rpt,phi,Ht,zot,bt)
    g3[i] = G_spiral_pipe(Fo,rpt,phin,Ht,zot,bt)
#    I = quad(g_integrand, phi1, phi2, args=(Fo,rpt,Ht,zot,bt),full_output = 1)
#    g2[i] = I[0]/(phi2-phi1)
gv2 = np.loadtxt('..\\data\\mean.txt')
g2 = gv2[:,1]
p1 = plot(Fov,g1,color = 'k',linewidth = 2,marker = 'o',markersize=8,label = 'z = H/2')
p2 = plot(Fov,g2,color = 'k',linewidth = 2,marker = 'x',markersize=8,label = 'mean')
p3 = plot(Fov,g3,color = 'k',linewidth = 2,marker = '+',markersize=8,label = 'z = H/4')
legend(fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14
ylabel(' $ k\Delta T/q\'$',fontname='Times new Roman',fontsize = fts)
xlabel(' Fo',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()