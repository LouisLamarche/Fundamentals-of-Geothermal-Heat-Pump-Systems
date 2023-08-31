#coding: utf-8
#
# Exemple 5.5 dimensionnemnt approche suédoise
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from time import *
# fichier de fonction g (Eskilson) tabuléees pour champ 2 x 2 pour b = 0.05,0.1,0.2,0.4,0.8
zo = 0
ro = 0.25/2
b = 0.5
cas = 'b'
if cas == 'a':
    H = 5
    tit = 'H = 5 m'
else:
    H = 50
    tit = 'H = 50 m'
Ht = H/ro
bt = b/ro
Fo = 5
rp = 0.01
rpt = rp/ro
zot = zo/ro
ome = 2*pi/b
omet = ome*ro
zv = np.linspace(0,H,30)
nz = len(zv)
th = np.zeros(nz)
for i in range(0,nz):
    z = zv[i]
    phi = ome*z
    th2 = G_spiral_pipe(Fo,rpt,phi,Ht,zot,bt)
    th[i] = th2

Fo = 500
thn = np.zeros(nz)
for i in range(0,nz):
    z = zv[i]
    phi = ome*z
    th2 = G_spiral_pipe(Fo,rpt,phi,Ht,zot,bt)
    thn[i] = th2
x = zv
p1 = plot(th,x,color = 'k',linewidth = 2,marker = 'o',markersize=8,label = 'Fo = 5')
p2 = plot(thn,x,color = 'k',linewidth = 2,marker = 'x',markersize=8,label = 'Fo = 500')
legend(fontsize = 14)
axis([0,1,H,0])
ax = gca()

grid(True,which='both')
fts = 16
ftst = 14
ylabel(' depth(m)',fontname='Times new Roman',fontsize = fts)
xlabel(' $ k\Delta T/q\'$',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
title(tit,fontsize = 22,fontname='Times new Roman')
tight_layout()
show()