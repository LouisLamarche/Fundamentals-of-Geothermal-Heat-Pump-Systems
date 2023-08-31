#coding: latin-1
#
# Exemple 5.5 dimensionnemnt approche suédoise
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from time import *
# fichier de fonction g (Eskilson) tabuléees pour champ 2 x 2 pour b = 0.05,0.1,0.2,0.4,0.8
zo = 0
rb = 0.25/2
b = 0.05
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
Fov = alm*t/rb**2
r = 0.25
rt = r/rb
zot = zo/rb
zt = z/rb
phi = 0
q = 400
qp = q/H
nt = len(Fov)
g1 = np.zeros(nt)
g2  = np.zeros(nt)
g3  = np.zeros(nt)
g4  = np.zeros(nt)
g5  = np.zeros(nt)
tic = time()
for i in range(0,nt):
    Fo = Fov[i]
    g1[i] = qp*G_function_fsc(Fo,rt,zt,Ht,zot)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nt):
    Fo = Fov[i]
    g2[i] =  qp*G_function_ring(Fo,rt,zt,Ht,zot,bt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nt):
    Fo = Fov[i]
    g3[i] =  qp*G_function_spiral(Fo,rt,zt,phi,Ht,zot,bt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nt):
    Fo = Fov[i]
    g4[i] =  qp*G_function_ils(Fo,rt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
tic = time()
for i in range(0,nt):
    Fo = Fov[i]
    g5[i] =  qp*G_function_ics(Fo,rt)/k
toc = time() - tic
print ('fin du calcul de la fonction g',toc )
x = np.log(Fov)
#x = t
p1 = plot(x,g1,color = 'k',linewidth = 2,label = 'Man''s solid cylinder')
p2 = plot(x,g2,color = 'k',linestyle = 'None',marker = 'o',markersize=8,label = 'Cui ring''s model')
p3 = plot(x,g3,color = 'k',linestyle = 'None',marker = 'x',markersize=8,label = 'Park''s spiral model')
p4 = plot(x,g4,color = 'k',linestyle = 'None',marker = '+',markersize=8,label = 'ILS')
p5 = plot(x,g5,color = 'k',linestyle = 'None',marker = 's',markersize=8,label = 'ICS')
legend(fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14

xlabel(' log(Fo)',fontname='Times new Roman',fontsize = fts)
ylabel(' $ \Delta T$',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()
