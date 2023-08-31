#coding: utf-8
#
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from time import *
# fichier de fonction g (Eskilson) tabuléees pour champ 2 x 2 pour b = 0.05,0.1,0.2,0.4,0.8
zo = 0
rb = 0.25/2
b = 0.5
H = 4
z = 1.8
Ht = H/rb
bt = b/rb
t = np.arange(0,1500,50)
Cp = 800
k = 0.26
rho = 1397
alp = k/(rho*Cp)
alm = alp*60
Fo = 1.5
rt = 1.4
zot = zo/rb
dh = H/60
zv = np.arange(0,H+2*dh,dh)
phi = 0
q = 400
qp = q/H
nz = len(zv)
g1 = np.zeros(nz)
g2  = np.zeros(nz)
g3  = np.zeros(nz)
tic = time()
for i in range(0,nz):
    zt = zv[i]/rb
    g1[i] = qp*G_function_fsc(Fo,rt,zt,Ht,zot)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nz):
    zt = zv[i]/rb
    g2[i] =  qp*G_function_ring(Fo,rt,zt,Ht,zot,bt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nz):
    zt = zv[i]/rb
    g3[i] =  qp*G_function_spiral(Fo,rt,zt,phi,Ht,zot,bt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
x = zv
#x = t
p1 = plot(g1,x,color = 'k',linewidth = 2,label = 'Man''s solid cylinder')
p2 = plot(g2,x,color = 'k',linewidth = 2,marker = 'o',markersize=8,label = 'Cui ring''s model')
p3 = plot(g3,x,color = 'k',linewidth = 2,marker = 'x',markersize=8,label = 'Park''s spiral model')
#p2 = plot(g2,x,color = 'k',linestyle = 'none',marker = 'o',markersize=8,label = 'Cui ring''s model')
#p3 = plot(g3,x,color = 'k',linestyle = 'none',marker = 'x',markersize=8,label = 'Park''s spiral model')
legend(fontsize = 14)
axis([0,80,5,0])
ax = gca()

grid(True,which='both')
fts = 16
ftst = 14

ylabel(' depth(m)',fontname='Times new Roman',fontsize = fts)
xlabel(' $ \Delta T$',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()
