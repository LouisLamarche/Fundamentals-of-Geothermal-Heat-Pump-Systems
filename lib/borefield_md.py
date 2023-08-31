from utilities_md import *
from geothermal_md import *
import numpy as np
from copy import deepcopy
import ctypes
from os import getcwd
from numpy.ctypeslib import ndpointer
current_dir = getcwd()
mes_fcts = ctypes.cdll.LoadLibrary(current_dir  + '\\my_fct2.so')
mon_prod  = mes_fcts.dfun
mon_prod.restype = None
mon_prod.argtypes = [ctypes.c_size_t,ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_size_t,ctypes.c_size_t,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_size_t,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

def Initialisation_Aplus(dtt,rr,nt,zb,T1,T2,nzi):


    ntt = 50
    nb = len(zb)
    D = np.zeros((nb,2))
    ntot = int((nb**2 + nb)/2)
    gam = 9*rr**2
    c = 1.8*np.sqrt(gam)
    zf = 15/np.sqrt(gam)
    uf = 1-1/(1+c*zf)
    du = uf/nzi
    u = np.arange(0,uf+du,du)
    zzz = (1/(1-u)-1)/c
    nz = len(zzz)-1
    y2 = np.zeros(nz)
    Dz = np.zeros(nz)
    uij = np.zeros((nz,ntot))
    xij = np.zeros((nb,nb))
    Dz = np.diff(zzz)
    zv = zzz[1:nz+1]
    Nga = 6
    t1b = np.log(2)/max(zv)**2
    t2b = np.log(2)/min(zv)**2*Nga
    t1a = 0.75*dtt
    t2a = 1.25*nt*dtt
    t1 = min(t1a,t1b)
    t2 = max(t2a,t2b)
    ttt = np.logspace(np.log10(t1),np.log10(t2),ntt)
    xxx = np.log(ttt)
    # calcul de uii
    Gfct = np.vectorize(g_function_fls)
    gg = Gfct(ttt,rbb = rr,zob=0)
    ders = der_sec(xxx,gg)
    inv_laplace = cls_inv_laplace(Nga)
    for i in range(0,nz):
        xq = zv[i]*zv[i]
        y2[i]=-2*zv[i]*inv_laplace.gavsteh_log_spline(xq,Nga,i,ders,xxx,gg)/(2*pi)
    ds = np.exp(-zv**2*dtt)
    dss = (1-ds)*y2
    for i in range(0,nb):
        k = int(nb*i - i*(i+1)/2 + i)
        uij[:,k] = dss
        xij[i,i] = np.sum(dss*Dz)
    for i in range(0,nb):
        for j in range(i+1,nb):
            rrr = np.linalg.norm(zb[i]-zb[j])
            gh =  Gfct(ttt,rbb = rrr,zob=0)
            ders = der_sec(xxx,gh)
            for  iz in range(0,nz):
                xq = (zv[iz])**2
                y2[iz]=-2*zv[iz]*inv_laplace.gavsteh_log_spline(xq,Nga,i,ders,xxx,gh)/(2*pi)
            dss = (1-ds)*y2
            k = int(nb*(i) - i*(i+1)/2 + j)
            uij[:,k] =  dss
            xij[i,j] = np.sum(dss*Dz)
            xij[j,i] = xij[i,j]
        for j in range(0,nb):
            if T1[j]> 0:
                D[i,0] = D[i,0] - xij[i,j]
            if T2[j]> 0:
                D[i,1] = D[i,1] - xij[i,j]
    return Dz,ds,uij,xij,D


def Initialisation_Aplus_new(dtt,rr,nt,zb,T1,T2,nz):


    nb = len(zb)
    Nga = 6
    xv,wv = lgwt(nz,0,1)
    D = np.zeros((nb,2))
    ntot = int((nb**2 + nb)/2)
    gam = 9*rr**2
    sgam = 3*rr
    wv = wv/xv**2/sgam
    yy  = (1/xv - 1)/sgam
    y2 = np.zeros(nz)
    uij = np.zeros((nz,ntot))
    xij = np.zeros((nb,nb))
    inv_laplace = cls_inv_laplace(Nga)
    for i in range(0,nz):
        xq = yy[i]*yy[i]
        y2[i]=-2*yy[i]*inv_laplace.gavsteh(xq,Nga,i,g_function_fls,rbb = rr)/(2*pi)
    ds = np.exp(-yy**2*dtt)
    uu = y2
    dss = (1-ds)*uu
    for i in range(0,nb):
        k = int(nb*i - i*(i+1)/2 + i)
        uij[:,k] = dss
        xij[i,i] = np.sum(dss*wv)
    for i in range(0,nb):
        for j in range(i+1,nb):
            rrr = np.linalg.norm(zb[i]-zb[j])
            for  iz in range(0,nz):
                xq = (yy[iz])**2
                y2[iz] = -2*yy[iz]*inv_laplace.gavsteh(xq,Nga,i,g_function_fls,rbb = rrr)/(2*pi)
            uu = y2
            dss = (1-ds)*uu
            k = int(nb*(i) - i*(i+1)/2 + j)
            uij[:,k] =  dss
            xij[i,j] = np.sum(dss*wv)
            xij[j,i] = xij[i,j]
        for j in range(0,nb):
            if T1[j]> 0:
                D[i,0] = D[i,0] - xij[i,j]
            if T2[j]> 0:
                D[i,1] = D[i,1] - xij[i,j]
    return wv,ds,uij,xij,D

def Initialisation_Aplus_new1(dtt,rr,nt,zb,nz):


    nb = len(zb)
    Nga = 6
    xv,wv = lgwt(nz,0,1)
    D = np.zeros((nb,1))
    ntot = int((nb**2 + nb)/2)
    gam = 9*rr**2
    sgam = 3*rr
#    wv = wv/xv**2/sgam
    yy  = (1/xv - 1)/sgam
    y2 = np.zeros(nz)
    uij = np.zeros((nz,ntot))
    xij = np.zeros((nb,nb))
    inv_laplace = cls_inv_laplace(Nga)
    for i in range(0,nz):
        xq = yy[i]*yy[i]
        y2[i]=-2*yy[i]*inv_laplace.gavsteh(xq,Nga,i,g_function_fls,rbb = rr)/(2*pi)
    ds = np.exp(-yy**2*dtt)
    uu = y2/xv**2/sgam
    dss = (1-ds)*uu
    for i in range(0,nb):
        k = int(nb*i - i*(i+1)/2 + i)
        uij[:,k] = dss
        xij[i,i] = np.sum(dss*wv)
    for i in range(0,nb):
        for j in range(i+1,nb):
            rrr = np.linalg.norm(zb[i]-zb[j])
            for  iz in range(0,nz):
                xq = (yy[iz])**2
                y2[iz] = -2*yy[iz]*inv_laplace.gavsteh(xq,Nga,i,g_function_fls,rbb = rrr)/(2*pi)
            uu = y2/xv**2/sgam
            dss = (1-ds)*uu
            k = int(nb*(i) - i*(i+1)/2 + j)
            uij[:,k] =  dss
            xij[i,j] = np.sum(dss*wv)
            xij[j,i] = xij[i,j]
        for j in range(0,nb):
            D[i,0] = D[i,0] - xij[i,j]
    return wv,ds,uij,xij,D


def Initialisation_Aplus_newb(dtt,rr,nt,zb,T1,T2,nz):


    nb = len(zb)
    Nga = 6
    xv,wv = lgwt(nz,0,1)
    D = np.zeros((nb,2))
    ntot = int((nb**2 + nb)/2)
    gam = 9*rr**2
    sgam = 3*rr
    xi = 1.25/sgam
    xf = 45/sgam
    xv1,wv1 = lgwt(nz,0,xi)
    xv2,wv2 = lgwt(nz,xi,xf)
    inv_laplace1 = cls_inv_laplace(Nga)
    inv_laplace2 = cls_inv_laplace(Nga)
    y1 = np.zeros(nz)
    y2 = np.zeros(nz)
    uij1 = np.zeros((nz,ntot))
    uij2 = np.zeros((nz,ntot))
    xij = np.zeros((nb,nb))
    for i in range(0,nz):
        xq1 = xv1[i]*xv1[i]
        xq2 = xv2[i]*xv2[i]
        y1[i]=-2*xv1[i]*inv_laplace1.gavsteh(xq1,Nga,i,g_function_fls,rbb = rr)/(2*pi)
        y2[i]=-2*xv2[i]*inv_laplace2.gavsteh(xq2,Nga,i,g_function_fls,rbb = rr)/(2*pi)
    ds1 = np.exp(-xv1**2*dtt)
    u1 = y1
    dss1 = (1-ds1)*u1
    ds2 = np.exp(-xv2**2*dtt)
    u2 = y2
    dss2 = (1-ds2)*u2
    for i in range(0,nb):
        k = int(nb*i - i*(i+1)/2 + i)
        uij1[:,k] = dss1
        uij2[:,k] = dss2
        xij[i,i] = np.sum(dss1*wv1)+np.sum(dss2*wv2)
    for i in range(0,nb):
        for j in range(i+1,nb):
            rrr = np.linalg.norm(zb[i]-zb[j])
            inv_laplace1 = cls_inv_laplace(Nga)
            inv_laplace2 = cls_inv_laplace(Nga)
            for  iz in range(0,nz):
                xq1 = xv1[iz]*xv1[iz]
                xq2 = xv2[iz]*xv2[iz]
                y1[iz]=-2*xv1[iz]*inv_laplace1.gavsteh(xq1,Nga,iz,g_function_fls,rbb = rrr)/(2*pi)
                y2[iz]=-2*xv2[iz]*inv_laplace2.gavsteh(xq2,Nga,iz,g_function_fls,rbb = rrr)/(2*pi)
            u1 = y1
            dss1 = (1-ds1)*u1
            ds2 = np.exp(-xv2**2*dtt)
            u2 = y2
            dss2 = (1-ds2)*u2
            k = int(nb*(i) - i*(i+1)/2 + j)
            uij1[:,k] =  dss1
            uij2[:,k] =  dss2
            xij[i,j] = np.sum(dss1*wv1) +  np.sum(dss2*wv2)
            xij[j,i] = xij[i,j]
        for j in range(0,nb):
            if T1[j]> 0:
                D[i,0] = D[i,0] - xij[i,j]
            if T2[j]> 0:
                D[i,1] = D[i,1] - xij[i,j]
    return wv1,ds1,uij1,wv2,ds2,uij2,xij,D



def calcul_Tb_Tf_plus(nb,Dz,nz,uij,thepp,CCL,A,Bn,ks,T1,T2,Fai):

    Fi = np.zeros((nz,nb))
    ntot = int((nb**2 + nb)/2)
    B = deepcopy(Bn)
    BB = sum(Dz[:,None]*Fai)/ks
    B[0:nb] = BB
    C = A @ B
    Tbi = C[0:nb]
    Tf1 = C[nb]
    Tf2 = C[nb+1]
    Tf = T1*Tf1 + T2*Tf2
    Tfout = thepp*Tf+(1-thepp)*Tbi
    qpi = CCL*(Tf-Tfout)
    sz = np.zeros(nz)
    for i in range(0,nb):
        mon_prod(i,uij,nz,ntot,qpi,nb,sz)
        Fi[:,i]  = Fai[:,i] + sz
    return Tfout,Fi,qpi


def calcul_Tb_Tf_plus_new(nb,wv,nz,uij,thepp,CCL,A,Bn,ks,T1,T2,Fai):

    Fi = np.zeros((nz,nb))
    ntot = int((nb**2 + nb)/2)
#    B = deepcopy(Bn)
    B = Bn.copy()
    BB = np.sum(wv[:,None]*Fai,axis = 0)/ks
    B[0:nb] = BB
    C = A @ B
    Tbi = C[0:nb]
    Tf1 = C[nb]
    Tf2 = C[nb+1]
    Tf = T1*Tf1 + T2*Tf2
    Tfout = thepp*Tf+(1-thepp)*Tbi
    qpi = CCL*(Tf-Tfout)
    sz = np.zeros(nz)
    for i in range(0,nb):
        mon_prod(i,uij,nz,ntot,qpi,nb,sz)
        Fi[:,i]  = Fai[:,i] + sz
    return Tfout,Fi,qpi,Tf1,Tf2

def calcul_Tb_Tf_plus_new2UT(nb,wv,nz,uij,XQ1,XQ2,A,Bn,ks,Fai):

    Fi = np.zeros((nz,nb))
    ntot = int((nb**2 + nb)/2)
    B = deepcopy(Bn)
    BB = sum(wv[:,None]*Fai)/ks
    B[0:nb] = BB
    C = A @ B
    Tbi = C[0:nb]
    Tbar = C[nb]
    Ttil = C[nb+1]
    qpi = (Ttil*XQ2+(Tbar-Tbi)*XQ1)
    sz = np.zeros(nz)
    for i in range(0,nb):
        mon_prod(i,uij,nz,ntot,qpi,nb,sz)
        Fi[:,i]  = Fai[:,i] + sz
    return Tbi,Fi,qpi,Tbar,Ttil


def calcul_Tb_Tf_plus_new1(nb,wv,nz,uij,thepp,CCL,A,Bn,ks,Fai):

    Fi = np.zeros((nz,nb))
    ntot = int((nb**2 + nb)/2)
    B = deepcopy(Bn)
    BB = sum(wv[:,None]*Fai)/ks
    B[0:nb] = BB
    C = A @ B
    Tbi = C[0:nb]
    Tf1 = C[nb]
    Tfout = thepp*Tf1+(1-thepp)*Tbi
    qpi = CCL*(Tf1-Tfout)
    sz = np.zeros(nz)
    for i in range(0,nb):
        mon_prod(i,uij,nz,ntot,qpi,nb,sz)
        Fi[:,i]  = Fai[:,i] + sz
    return Tfout,Fi,qpi

def calcul_Tb_Tf_plus_newb(nb,wv1,wv2,nz,uij1,uij2,thepp,CCL,A,Bn,ks,T1,T2,Fai1,Fai2):

    Fi1 = np.zeros((nz,nb))
    Fi2 = np.zeros((nz,nb))
    ntot = int((nb**2 + nb)/2)
    B = deepcopy(Bn)
    BB = sum(wv1[:,None]*Fai1)/ks + sum(wv2[:,None]*Fai2)/ks
    B[0:nb] = BB
    C = A @ B
    Tbi = C[0:nb]
    Tf1 = C[nb]
    Tf2 = C[nb+1]
    Tf = T1*Tf1 + T2*Tf2
    Tfout = thepp*Tf+(1-thepp)*Tbi
    qpi = CCL*(Tf-Tfout)
    sz1 = np.zeros(nz)
    sz2 = np.zeros(nz)
    for i in range(0,nb):
        mon_prod(i,uij1,nz,ntot,qpi,nb,sz1)
        Fi1[:,i]  = Fai1[:,i] + sz1
        mon_prod(i,uij2,nz,ntot,qpi,nb,sz2)
        Fi2[:,i]  = Fai2[:,i] + sz2
    return Tfout,Fi1,Fi2,qpi


def  Calcul_Q1_Q4(Rd1,Rd12,Rd13,Rd14,CC1,CC2):

    R1s = CC1*Rd1
    R12s = CC1*Rd12
    R13s = CC1*Rd13
    R14s = CC1*Rd14
    a = 1/R1s + 1/R12s + 1/R13s+ 1/R14s
    b = -1/R12s
    c = -1/R13s
    d = -1/R14s
    al = CC1/CC2
    h1 = a**2+ b**2 + c**2 - d**2
    h2 = a**2+ b**2 - c**2 + d**2
    h3 = a**2 - b**2 + c**2 + d**2
    h4 = -a**2 + b**2 + c**2 + d**2
    f1 = a + c - b - d
    al2 = al*al
    X = (a**2-c**2)*(al2+1)+2*al*(b**2-d**2)
    Y = np.sqrt(((a**2-c**2)*(al2+1)+2*al*(b**2-d**2))**2-4*al2*((a**2-b**2)**2+(c**2-d**2)**2-2*(a**2+b**2)*(c**2+d**2) + 8*a*b*c*d))
    gam2 = (X+Y)/2
    eta2 = (X-Y)/2
    gam = np.sqrt(gam2)
    eta = np.sqrt(eta2)
    den = gam2 - eta2
    den1 = den*gam
    den2 = den*eta
    G12i= gam2 + al2*(c**2-a**2) + al*(d**2-b**2) +(c*d-a*b)*(1+al)
    G32i =  (c*b-a*d)*(1-al)
    G14i= -1*(eta2 + al2*(c**2-a**2) + al*(d**2-b**2)) + (-c*d+a*b)*(1+al)
    G34i =  - (c*b-a*d)*(1-al)
    G42i =  (a*d - c*b)*(al-1)*al
    G22i = -gam2 + (a**2-c**2) + al2*(a*b-c*d)+ al*(a*b-c*d + b**2 - d**2)
    G24i = eta2 - (a**2-c**2) - al2*(a*b-c*d)- al*(a*b-c*d + b**2 - d**2)
    G44i =   - (a*d - c*b)*(al-1)*al
    G12p=0
    G12pp=(b*c-a*d)*(al-1)
    G12pp = G12pp/den
    G13p= -c*gam2+c*al2*h2 - 2*al2*a*b*d
    G13p = G13p/den1
    G13pp= -d*gam2 + al*d*h1 - 2*al*a*b*c
    G13pp = G13pp/den1
    G14p=0
    G14pp=-(b*c-a*d)*(al-1)
    G14pp = G14pp/den
    G15p= c*eta2 -c*al2*h2 + 2*al2*a*b*d
    G15p = G15p/den2
    G15pp= d*eta2 - al*d*h1 + 2*al*a*b*c
    G15pp = G15pp/den2
    G22p = -al*(b*c-a*d)*(al-1)
    G22p = G22p/den
    G22pp = 0
    G23p=  -al*d*gam2 + al2*(d*h1 - 2*a*b*c)
    G23p = G23p/den1
    G23pp=  -al*c*gam2 + al*(c*h2 - 2*a*b*d)
    G23pp = G23pp/den1
    G24p = al*(b*c-a*d)*(al-1)
    G24p = G24p/den
    G24pp = 0
    G25p=  al*d*eta2 - al2*(d*h1 - 2*a*b*c)
    G25p = G25p/den2
    G25pp=  al*c*eta2 - al*(c*h2 - 2*a*b*d)
    G25pp = G25pp/den2
    G32p=  gam2 - al2*(a**2-c**2) -al*(b**2-d**2)
    G32p = G32p/den
    G32pp=  (a*b-c*d)*(1+al)
    G32pp = G32pp/den
    G33p=  a*gam2 + al2*(a*h4-2*b*c*d)
    G33p = G33p/den1
    G33pp = b*gam2 + al*(b*h3-2*a*c*d)
    G33pp = G33pp/den1
    G34p=  -eta2 + al2*(a**2-c**2) + al*(b**2-d**2)
    G34p = G34p/den
    G34pp=  -(a*b-c*d)*(1+al)
    G34pp = G34pp/den
    G35p=  -a*eta2 - al2*(a*h4-2*b*c*d)
    G35p = G35p/den2
    G35pp = -b*eta2 - al*(b*h3-2*a*c*d)
    G35pp = G35pp/den2
    G42p=  al*(a*b-c*d)*(1+al)
    G42p = G42p/den
    G42pp=  gam2 + al*(d**2-b**2) + c**2 - a**2
    G42pp = G42pp/den
    G43p=  al*b*gam2 + al2*(b*h3-2*a*c*d)
    G43p = G43p/den1
    G43pp = al*a*gam2 + al*(a*h4-2*b*c*d)
    G43pp = G43pp/den1
    G44p=  -al*(a*b-c*d)*(1+al)
    G44p = G44p/den
    G44pp=  -eta2 - al*(d**2-b**2) + a**2 - c**2
    G44pp = G44pp/den
    G45p=  -al*b*eta2 - al2*(b*h3-2*a*c*d)
    G45p = G45p/den2
    G45pp = -al*a*eta2 - al*(a*h4-2*b*c*d)
    G45pp = G45pp/den2
    A11=(G12p-G32p)*np.cosh(gam)+(G13p-G33p)*np.sinh(gam)+(G14p-G34p)*np.cosh(eta)+(G15p-G35p)*np.sinh(eta)
    A12=(G12pp-G32pp)*np.cosh(gam)+(G13pp-G33pp)*np.sinh(gam)+(G14pp-G34pp)*np.cosh(eta)+(G15pp-G35pp)*np.sinh(eta)
    A21=(G22p-G42p)*np.cosh(gam)+(G23p-G43p)*np.sinh(gam)+(G24p-G44p)*np.cosh(eta)+(G25p-G45p)*np.sinh(eta)
    A22=(G22pp-G42pp)*np.cosh(gam)+(G23pp-G43pp)*np.sinh(gam)+(G24pp-G44pp)*np.cosh(eta)+(G25pp-G45pp)*np.sinh(eta)
    delta = (A11*A22 - A12*A21)
    I1 = al2*(c*h2 - a*h4 + 2*(b*c*d - a*b*d)) + al*(b*h3 - d*h1 + 2*(a*b*c-a*c*d)) + gam2*(b-a-c+d)
    I2 = -al2*(c*h2 - a*h4 + 2*(b*c*d - a*b*d)) - al*(b*h3 - d*h1 + 2*(a*b*c-a*c*d)) - eta2*(b-a-c+d)
    I3 = al2*(d*h1 - b*h3 + 2*(a*c*d - a*b*c)) + al*(a*h4 - c*h2 + 2*(a*b*d-b*c*d)) - al*gam2*(b-a-c+d)
    I4 = -al2*(d*h1 - b*h3 + 2*(a*c*d - a*b*c)) - al*(a*h4 - c*h2 + 2*(a*b*d-b*c*d)) + al*eta2*(b-a-c+d)
    I5 = 2/R1s*(f1*al*(al*(a-c) - b + d) - gam2)
    I6 = -2/R1s*(f1*al*(al*(a-c) - b + d) - eta2)
    I7 = -2*al/R1s*(f1*(al*(b-d) - a + c) + gam2)
    I8 = 2*al/R1s*(f1*(al*(b-d) - a + c) + eta2)
    Q1 = (A22*(I5*np.sinh(gam)/gam + I6*np.sinh(eta)/eta)- A12*(I7*np.sinh(gam)/gam + I8*np.sinh(eta)/eta))/(den*delta)
    Q2 = (A22*((G32i-G12i)*np.cosh(gam) + (G34i-G14i)*np.cosh(eta) - I1*np.sinh(gam)/gam - I2 *np.sinh(eta)/eta) - \
          A12*((G42i-G22i)*np.cosh(gam) + (G44i-G24i)*np.cosh(eta) - I3*np.sinh(gam)/gam - I4 *np.sinh(eta)/eta))/(den*delta)
    Q3 = (-A21*(I5*np.sinh(gam)/gam + I6*np.sinh(eta)/eta) + A11*(I7*np.sinh(gam)/gam + I8*np.sinh(eta)/eta))/(den*delta)
    Q4 = (-A21*((G32i-G12i)*np.cosh(gam) + (G34i-G14i)*np.cosh(eta) - I1*np.sinh(gam)/gam - I2 *np.sinh(eta)/eta) + \
        A11*((G42i-G22i)*np.cosh(gam) + (G44i-G24i)*np.cosh(eta) - I3*np.sinh(gam)/gam - I4 *np.sinh(eta)/eta))/(den*delta)
    return Q1,Q2,Q3,Q4
