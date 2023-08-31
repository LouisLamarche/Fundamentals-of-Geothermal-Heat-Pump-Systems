#coding: utf-8
#
#
from geothermal_md import *
import numpy as np
from numpy.fft import fft,rfft,ifft
from  matplotlib.pyplot import *
data = np.loadtxt("..\\data\\g0401n.txt")  #
# fichier de fonction g (Eskilson) tabuléees pour champ 2 x 2 pour b = 0.05,0.1,0.2,0.4,0.8
xv = data[:,0]  # x = log(t/ts)
bv = data[:,1]  # bv = 0.05,0.1,0.2,0.4,0.8
gv = data[:,2]  # fonction g pour les 5 valeurs de b
nw = len(gv)
ne = 5
bb = bv[0:ne]
# premier cas
i = 0   # b = 0.05
b = bb[i]
x_es = xv[0:nw:ne]
g_es = gv[i:nw:ne]       # fonction g tabulées de Eskilson pour le premier espacement

#
# generaton de la fonction g analytique par la LSF
nx = 4
ny = 1
nb = nx*ny
ns = 2
ntot = nb*ns

#
#
zt = np.zeros([nb,2])    # zt represente une matrice ayant les coordonnées adimensionnels ( b = B/H) des nb puits
k = 0
for i1  in range(0,nx):
    dx = i1*b
    for i2 in range(0,ny):
        dy = i2*b
        zt[k] = [dx,dy]
        k = k+1
#
#
#
nn = len(x_es)
rr = 0.0005
ala = 1e-6*3600*24*365
tc = 100**2/(9*ala)
dy = 12*4
nt = 20
ty = np.arange(0,nt*dy+dy,dy)
tt = ty/12/tc
tt = np.array([0,.05,.1,.15,.2,.25,.3])
tmax = max(tt)
nt = len(tt)
n1 = nt//2 + 1
n2 = (nt-1)//2
sig = 2*np.log(nt)/tmax
gg1 = np.zeros(nt)
for i in range (0,nt):
    gg1[i] = compute_g_function(zt,tt[i],rbb = rr)    # SI on n'ecrit pas D/H , la fonction prend 0.04 par défaut ( Comme Eskilson)
#
#
# Exemple d'utilisation de la fonction de Massimo Cimmino
# References
#    ----------
#    .. [#CimminoBernier2014] Cimmino, M., & Bernier, M. (2014). A
#       semi-analytical method to generate g-functions for geothermal bore
#       fields. International Journal of Heat and Mass Transfer, 70, 641-650.

nsy = 4
pr = 0
go = np.zeros(nt)
g1 = np.zeros(nt)
g2 = np.zeros(nt)
g3 = np.zeros(nt)
g4 = np.zeros(nt)
g5 = np.zeros(nt)
g6 = np.zeros(nt)
g7 = np.zeros(nt)
g8 = np.zeros(nt)
g9 = np.zeros(nt)
g10 = np.zeros(nt)
g11 = np.zeros(nt)
g12 = np.zeros(nt)
g13 = np.zeros(nt)
g14 = np.zeros(nt)
g15 = np.zeros(nt)
rr = 0.0005
z1  = 0.04
b = 0.05
h1 = 0.5
h2 = 0.5
z2 = z1 + h1
for i in range(0,nt):
    go[i] =  fct_hij(tt[i],rr,z1,z1,h1,h1)
    g1[i] =  fct_hij(tt[i],b,z1,z1,h1,h1)
    g2[i] =  fct_hij(tt[i],2*b,z1,z1,h1,h1)
    g3[i] =  fct_hij(tt[i],3*b,z1,z1,h1,h1)
    g4[i] =  fct_hij(tt[i],rr,z2,z1,h1,h1)
    g5[i] =  fct_hij(tt[i],b,z2,z1,h1,h1)
    g6[i] =  fct_hij(tt[i],2*b,z2,z1,h1,h1)
    g7[i] =  fct_hij(tt[i],3*b,z2,z1,h1,h1)
    g8[i] =  fct_hij(tt[i],rr,z1,z2,h1,h1)
    g9[i] =  fct_hij(tt[i],b,z1,z2,h1,h1)
    g10[i] =  fct_hij(tt[i],2*b,z1,z2,h1,h1)
    g11[i] =  fct_hij(tt[i],3*b,z1,z2,h1,h1)
    g12[i] =  fct_hij(tt[i],rr,z2,z2,h1,h1)
    g13[i] =  fct_hij(tt[i],b,z2,z2,h1,h1)
    g14[i] =  fct_hij(tt[i],2*b,z2,z2,h1,h1)
    g15[i] =  fct_hij(tt[i],3*b,z2,z2,h1,h1)
gs = np.zeros((nt,nsy,nsy))
gs[:,0,0]  = go + g3
gs[:,0,1]  = g4 + g7
gs[:,0,2]  = g1 + g2
gs[:,0,3]  = g5 + g6
gs[:,1,0]  = g8 + g11
gs[:,1,1]  = g12 + g15
gs[:,1,2]  = g9 + g10
gs[:,1,3]  = g13 + g14
gs[:,2,0]  = g1 + g2
gs[:,2,1]  = g5 + g6
gs[:,2,2]  = go + g1
gs[:,2,3]  = g4 + g5
gs[:,3,0]  = g9 + g10
gs[:,3,1]  = g13 + g14
gs[:,3,2]  = g8 + g9
gs[:,3,3]  = g12 + g13
gbs = np.zeros((nt,nsy,nsy),dtype=complex)
qpbs = np.zeros((nsy,nt),dtype=complex)
AA = np.zeros((nsy+1,nsy+1),dtype=complex)
ggcs = np.zeros(nt,dtype=complex)


if pr:
    print(array2string(gs[:,0,0],precision = 4,separator = '&'))
    print(array2string(gs[:,0,1] ,precision = 4,separator = '&'))
    print(array2string(gs[:,1,1] ,precision = 4,separator = '&'))

#
# Calcul de T de Laplcae des facteurs de réponse thermique
#
tmax = max(tt)
sig = 2*np.log(nt)/tmax
for i in range(0,nsy):
    for j in range(0,nsy):
        gab = fft(gs[:,i,j]*np.exp(-sig*tt))
        gbs[:,i,j] = gab
if pr:
    print(array2string(gbs[:,0,0],precision = 4,separator = '&'))
    print(array2string(gbs[:,0,1] ,precision = 4,separator = '&'))
    print(array2string(gbs[:,1,1] ,precision = 4,separator = '&'))

B = np.zeros(nt)
B[0] = ntot
Bb =  fft(B*np.exp(-sig*tt))
BB = np.zeros(nsy+1,dtype=complex)
BB[nsy] = Bb[0]

for i in range(0,nt):
    A = gbs[i,:,:]
    AA[0:nsy,0:nsy] = A
    AA[nsy,0:nsy] = 2
    AA[0:nsy,nsy] = -1;
    HH = np.linalg.solve(AA, BB)
    ggcs[i] = HH[nsy]
    qpbs[:,i] = HH[0:nsy]
if pr:
    print(array2string(ggcs,precision = 4,separator = '&'))
    print(array2string(qpbs[0,:] ,precision = 4,separator = '&'))
    print(array2string(qpbs[1,:] ,precision = 4,separator = '&'))
ggc = np.real(np.exp(sig*tt)*ifft(ggcs))
qpc = np.zeros((nsy,nt-1))
for i in range(0,nsy):
   B = qpbs[i,:]
   Bi =  np.real(np.exp(sig*tt)*ifft(B))
   qpc[i,0] = Bi[0]
   for j in range(1,nt-1):
        qpc[i,j] =  qpc[i,j-1] + Bi[j]
ggc[0] = 0
if pr:
    print(array2string(ggc,precision = 4,separator = '&'))
    print(array2string(gg1,precision = 4,separator = '&'))

H = 100
B = b*H
r_b = rr*H
zo = 0.04*H
# LA fonction demande des valeurs dimensionnelles
#reField = rectangle_field(nx, ny, B,B,H,zo,r_b)
# Calculate the g-function for uniform borehole wall temperature
#nSegments = 2
#alpha = 1e-6
#ts = H**2/(9.*alpha)
ttt = tt[1:nt]
#time = ttt*ts

#gg1_Cimmino,qp_c = uniform_temperature(boreField, time, alpha,nSegments=nSegments,disp=False)
#

# plot

print('q11',np.array2string(qpc[0,:], formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))
print('q12',np.array2string(qpc[1,:], formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))
print('q21',np.array2string(qpc[2,:], formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))
print('q22',np.array2string(qpc[3,:], formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))

print('Type I',np.array2string(gg1, formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))
print('Type III',np.array2string(ggc, formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))

xx = np.log(ttt)
gges  = np.interp(xx,x_es,g_es)
print('Tabulated values',np.array2string(gges, formatter={'float_kind':lambda x: "%.4f" % x},separator = '\t'))

