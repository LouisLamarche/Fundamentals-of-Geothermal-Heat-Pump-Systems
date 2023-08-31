import numpy as np
from scipy.integrate import nquad,quad
from scipy.special import expn,j0,j1,y0,y1,k0,i0,i1
#from collections import *
from scipy.optimize import fsolve,minimize,newton,brentq
from scipy.interpolate import interpn
from scipy.special import factorial
from CoolProp.CoolProp import *
from numba import float64,jit
from math import erfc,erf,comb
import warnings
warnings.filterwarnings("ignore")

hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])              # hours in each month
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])   # cumulate hours

pi = np.pi
Inf = np.Inf



def calcul_kk(nr,jj,k):
    j = k + jj
    if j < nr:
       kk = j
    else:
       kk = 2*nr-(j+2)
    return kk




def g_function_spiral_app(Fo,xoi,yoi,xoj,yoj,rpt,h):
    dij = np.sqrt((xoi  - xoj)**2 + (yoi  - yoj)**2)
    dij_i = np.sqrt(dij**2 + 4*h*h)
    y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
    teta = pi*y
    return teta
def G_function_slinky_app(Fo,xoi,yoi,xoj,yoj,rpt,h):
    dij = np.sqrt((xoi  - xoj)**2 + (yoi  - yoj)**2)
    dij_i = np.sqrt(dij**2 + 4*h*h)
    y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
    teta = y/(2)
    return teta


@jit(float64(float64, float64, float64, float64, float64, float64, float64, float64, float64), nopython=True, cache=True)
def  g_integrand_xiong(ome,phi,Fo,xoi,yoi,xoj,yoj,rpt,h):
    d1 = np.sqrt((xoi + (1 - rpt)*np.cos(phi) - xoj - np.cos(ome))**2 + (yoi + (1 - rpt)*np.sin(phi) - yoj - np.sin(ome))**2)
    d2 = np.sqrt((xoi + (1 + rpt)*np.cos(phi) - xoj - np.cos(ome))**2 + (yoi + (1 + rpt)*np.sin(phi) - yoj - np.sin(ome))**2)
    dij = (d1 + d2)/2.0
    dij_i = np.sqrt(dij**2 + 4*h*h)
    y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
    return y


def g_function_spiral(Fo,xoi,yoi,xoj,yoj,rpt,h):
#
# fonction g Xiong
# rt  = r/ro
# ht= h/ro
# zt = z/ro

    zet = np.arctan2(yoj - yoi,xoj - xoi)
    I = nquad(g_integrand_xiong, [[0, 2*pi],[zet, zet+pi]],args =(Fo,xoi,yoi,xoj,yoj,rpt,h))
#    I = nquad(g_integrand_xiong, [[0, 2*pi],[0, 2*pi]],args =(Fo,xoi,yoi,xoj,yoj,rpt,h))
    teta = I[0]/(2*pi)
    return teta

def G_function_slinky_u(Fo,xoi,yoi,xoj,yoj,rpt,h):
#
# fonction g Xiong
# rt  = r/ro
# ht= h/ro
# zt = z/ro

    zet = np.arctan2(yoj - yoi,xoj - xoi)
    I = nquad(g_integrand_xiong, [[0, 2*pi],[zet, zet+pi]],args =(Fo,xoi,yoi,xoj,yoj,rpt,h))
#    I = nquad(g_integrand_xiong, [[0, 2*pi],[0, 2*pi]],args =(Fo,xoi,yoi,xoj,yoj,rpt,h))
    teta = I[0]/(4*pi**2)
    return teta

def calcul_fonction_g_xiong(Rb,Fo,rpt,zt,nt,ntr,r1=7,r2=25):
    nb = len(Rb)
    nsym = nt - 1
    dij = np.zeros(nsym)
    Qij1 = np.zeros(nsym)
    zi = Rb[0,:]
    xoi = zi[0]
    yoi = zi[1]
    for j in range(1,nt):
        zj = Rb[j,:]
        xoj = zj[0]
        yoj = zj[1]
        dij = np.linalg.norm(zi-zj)
        if dij > r2 + 2:
            Q = 0.0
        elif dij > r1 + 2:
            Q = g_function_spiral_app(Fo,xoi,yoi,xoj,yoj,rpt,zt)
        else:
            Q = g_function_spiral(Fo,xoi,yoi,xoj,yoj,rpt,zt)
        Qij1[j-1] = Q
    if ntr > 1:
        Qij2 = np.zeros(nt)
        zi = Rb[0,:]
        xoi = zi[0]
        yoi = zi[1]
        for j in range(0,nt):
            jj = j+nt
            zj = Rb[jj,:]
            xoj = zj[0]
            yoj = zj[1]
            dij = np.linalg.norm(zi-zj)
            if dij > r2 + 2:
                Q = 0.0
            elif dij > r1 + 2 :
                Q = g_function_spiral_app(Fo,xoi,yoi,xoj,yoj,rpt,zt)
            else:
                Q = g_function_spiral(Fo,xoi,yoi,xoj,yoj,rpt,zt)
            Qij2[j] = Q
    gb = np.zeros(nb)
    gbb = np.zeros((nb,nb))
    Qi = g_function_spiral(Fo,0,0,0,0,rpt,zt)
    for jj in range(0,nt):
        gbb[jj,jj] = Qi
        for k in range(jj+1,nt):
            Q =  Qij1[k-jj-1]
            gbb[jj,k]= Q
            gbb[k,jj] = Q
    if ntr > 1:
        for jj in range(0,nt):
            for k in range(0,nt):
                nj = nt+k
                kk =  calcul_kk(nt,jj,k)
                Q =  Qij2[kk]
                gbb[jj,nj ] = Q
                gbb[nj,jj] = Q
        for j in range(0,nt):
            jj = j+nt
            gbb[jj,jj] = Qi
            for k in range(j+1,nt):
                kk = k+nt
                Q =  Qij1[k-j-1]
                gbb[jj,kk]= Q
                gbb[kk,jj] = Q
    for jj in range(0,nb):
        gb[jj] = np.sum(gbb[jj,:])
    g2 = np.mean(gb)
    return g2

def G_function_slinky(Rb,Fo,rpt,zt,pt,ntr = 1,r1=7,r2=25):
    nb = len(Rb)
    nt = int(nb/ntr)
    nsym = nt - 1
    dij = np.zeros(nsym)
    Qij1 = np.zeros(nsym)
    zi = Rb[0,:]
    xoi = zi[0]
    yoi = zi[1]
    for j in range(1,nt):
        zj = Rb[j,:]
        xoj = zj[0]
        yoj = zj[1]
        dij = np.linalg.norm(zi-zj)
        if dij > r2 + 2:
            Q = 0.0
        elif dij > r1 + 2:
            Q = G_function_slinky_app(Fo,xoi,yoi,xoj,yoj,rpt,zt)
        else:
            Q = G_function_slinky_u(Fo,xoi,yoi,xoj,yoj,rpt,zt)
        Qij1[j-1] = Q
    if ntr > 1:
        Qij2 = np.zeros(nt)
        zi = Rb[0,:]
        xoi = zi[0]
        yoi = zi[1]
        for j in range(0,nt):
            jj = j+nt
            zj = Rb[jj,:]
            xoj = zj[0]
            yoj = zj[1]
            dij = np.linalg.norm(zi-zj)
            if dij > r2 + 2:
                Q = 0.0
            elif dij > r1 + 2 :
                Q = G_function_slinky_app(Fo,xoi,yoi,xoj,yoj,rpt,zt)
            else:
                Q = G_function_slinky_u(Fo,xoi,yoi,xoj,yoj,rpt,zt)
            Qij2[j] = Q
    gb = np.zeros(nb)
    gbb = np.zeros((nb,nb))
    Qi = G_function_slinky_u(Fo,0,0,0,0,rpt,zt)
    for jj in range(0,nt):
        gbb[jj,jj] = Qi
        for k in range(jj+1,nt):
            Q =  Qij1[k-jj-1]
            gbb[jj,k]= Q
            gbb[k,jj] = Q
    if ntr > 1:
        for jj in range(0,nt):
            for k in range(0,nt):
                nj = nt+k
                kk =  calcul_kk(nt,jj,k)
                Q =  Qij2[kk]
                gbb[jj,nj ] = Q
                gbb[nj,jj] = Q
        for j in range(0,nt):
            jj = j+nt
            gbb[jj,jj] = Qi
            for k in range(j+1,nt):
                kk = k+nt
                Q =  Qij1[k-j-1]
                gbb[jj,kk]= Q
                gbb[kk,jj] = Q
    for jj in range(0,nb):
        gb[jj] = np.sum(gbb[jj,:])
    g2 = pt*np.mean(gb)/(2*pi+2*pt)
    return g2


@jit(float64(float64), nopython=True, cache=True)
def ierf(x):
    y = x*erf(x) -(1-np.exp(-x*x))/np.sqrt(pi)
    return y


def g_integrand(x,Fo,rt):
    num = (np.exp(-Fo*x*x)-1)*(j0(rt*x)*y1(x)-y0(rt*x)*j1(x))
    den = x*x*(j1(x)*j1(x) + y1(x)*y1(x))
    y = num/den
    return y

def G_function_ics(Fo,rt=1):
    def g_integranda(x,Fo):
        num = (1-np.exp(-Fo*x*x))
        den = x*x*x*(j1(x)*j1(x) + y1(x)*y1(x))
        y = 2*num/den/pi
        return y
    if rt == 1:
        I = quad(g_integranda, 0, Inf, args=(Fo))
    else:
        I = quad(g_integrand, 0, Inf, args=(Fo,rt),full_output = 1)
    G = I[0]/(pi*pi)
    return G


def g_integrand_st(x,Fo,delt,gam,kt,nu,R):
    d1 = kt*gam
    z1 = x*gam
    z2  = x*(delt*gam)
    Rx = 1.0 - nu*R*x*x
    A = Rx*y1(z2) - x*nu*y0(z2)/(delt*d1)
    B = Rx*j1(z2) - x*nu*j0(z2)/(delt*d1)
    phi =  A*(y0(x)*j1(z1)*d1-y1(x)*j0(z1))+B*(y1(x)*y0(z1)-y0(x)*y1(z1)*d1)
    psi =  B*(j0(x)*y1(z1)*d1-j1(x)*y0(z1))+A*(j1(x)*j0(z1)-j0(x)*j1(z1)*d1)
    num = 1.0-np.exp(-Fo*x*x)
    den = (phi*phi+psi*psi)*x*x*x*x*x*(delt*delt*gam*gam)
    y = num/den
    return y

def G_function_st(Fo,ret,gam,kt,nu,R):
    I = quad(g_integrand_st, 0, Inf, args=(Fo,ret,gam,kt,nu,R),full_output = 1)
    G = I[0]*8/(pi*pi*pi*pi*pi)
    return G

def g_function_st(tt,rr,ret,gam,kt,nu,R):
    Fo = tt/(3*rr)**2
    req = ret/rr
    G = G_function_st(Fo,req,gam,kt,nu,R)
    g = G*2*pi
    return g

def G_function_isca(Fo,rt = 1):
    def g_integrand_mana(x,Fo):
        z = 1/(Fo -x)
        z2 = 0.5*rt*z
        z3 = (rt**2 + 1)/4*z
        if z2 < 20.0:
            y = z*np.exp(-z3)*i0(z2)
        else:
            y = np.sqrt(2*z2/pi)*np.exp(-z2*(rt-1)**2/(2*rt))*(1+1/(8*z2)+9/(128*z2**2))
        return y
    I = quad(g_integrand_mana, 0, Fo, args=(Fo),full_output = 1)
    G = I[0]/(4*pi)
    return G

def G_function_isc(Fo,rt = 1):
    def g_integrand_manb(x,Fo):
        z = (rt**2 + 1 - 2*rt*np.cos(x))/(4*Fo)
        y = expn(1,z)
        return y
    I = quad(g_integrand_manb, 0, pi, args=(Fo),full_output = 1)
    G = I[0]/(4*pi*pi)
    return G

def G_function_fsc(Fo,r,z,h,zo):
    def g_integrand_fsc(x,Fo,r,z,h,zo):

        de = 2*np.sqrt(Fo-x)
        a1 = (r**2+1)/(de*de)
        a2= 2*r/(de*de)
        beta1 =(z-zo)/de
        beta2 = (z+zo)/de
        beta3 = (z-h-zo)/de
        beta4 = (h+zo+z)/de
        y2 = erf(beta1)+erf(beta2)-(erf(beta3)+erf(beta4))
        y1 = np.exp(-a1)/(Fo-x)
        if a2 < 20.0:
            y3 = y1*i0(a2)
        else:
            y3 = np.sqrt(2*a2/pi)*np.exp(-a2*(r-1)**2/(2*r))*(1+1/(8*a2)+9/(128*a2**2))
        y = y2*y3
        return y
    I = quad(g_integrand_fsc, 0, Fo, args=(Fo,r,z,h,zo),full_output = 1)
    G = I[0]/(8*pi)
    return G



def G_function_ringc(Fo,r,z,h,zo,b):
    m = int(h/b)
    def g_integrand_ringc(x,Fo,r,z,h,zo,b):

        de = 2*sqrt(Fo-x)
        a1 = (r**2+1)/(de*de)
        a2= 2*r/(de*de)
        beta1 =(z-zo)/de
        beta2 = (z+zo)/de
        y2 = 0
        for i in range(0,m):
            zi = zo + (i + 0.5)*b
            beta1 = (z-zi)/de
            beta2 = (z+zi)/de
            y2 = y2 +  np.exp(-beta1**2) - np.exp(-beta2**2)
        y1 = np.exp(-a1)/(Fo-x)**1.5
        if a2 < 20.0:
            y3 = y1*i0(a2)
        else:
            y3 = 1/(np.sqrt(pi*r)*(Fo-x))*np.exp(-a2*(r-1)**2/(2*r))*(1+1/(8*a2)+9/(128*a2**2))
        y = y2*y3
        return y
    I = quad(g_integrand_ringc, 0, Fo, args=(Fo,r,z,h,zo,b),full_output = 1)
    G = b*I[0]/(8*pi*sqrt(pi))
    return G
def G_function_ring(Fo,r,z,h,zo,b):
    m = int(h/b)
    def g_integrand_ring(x,Fo,r,z,h,zo,b):

        d = r**2 + 1 - 2*r*np.cos(x)
        a1 = 2*np.sqrt(Fo)
        y2 = 0
        for i in range(0,m):
            zi = zo + (i + 0.5)*b
            b1 = np.sqrt(d + (z-zi)**2)
            b2 = np.sqrt(d + (z+zi)**2)
            y2 = y2 +  erfc(b1/a1)/b1 - erfc(b2/a1)/b2
        return y2
    I = quad(g_integrand_ring, 0, pi, args=(Fo,r,z,h,zo,b),full_output = 1)
    G = b*I[0]/(4*pi**2)
    return G

def G_function_spiral(Fo,r,z,phi,h,zo,b):
    om = 2*pi/b
    def g_integrand_spiral(x,Fo,r,z,phi,h):
        am = np.sqrt(r**2 + 1 - 2*r*np.cos(phi - om*x) +  (z - x)**2)
        ap = np.sqrt(r**2 + 1 - 2*r*np.cos(phi - om*x) +  (z + x)**2)
        y = erfc(am/(2*np.sqrt(Fo)))/am - erfc(ap/(2*np.sqrt(Fo)))/ap
        return y
    I = quad(g_integrand_spiral, zo, zo+h, args=(Fo,r,z,phi,h),full_output = 1)
    G = I[0]/(4*pi)
    return G


def G_spiral_pipe(Fo,rpt,phi,Ht,zot,bt):
    omet = 2*pi/bt
    theta = 0
    rr = 1 + rpt*np.cos(theta)
    zz = phi/omet + rpt*np.sin(theta)
    y1 = G_function_spiral(Fo,rr,zz,phi,Ht,zot,bt)
    theta = pi/2
    rr = 1 + rpt*np.cos(theta)
    zz = phi/omet + rpt*np.sin(theta)
    y2 = G_function_spiral(Fo,rr,zz,phi,Ht,zot,bt)
    theta = pi
    rr = 1 + rpt*np.cos(theta)
    zz = phi/omet + rpt*np.sin(theta)
    y3 = G_function_spiral(Fo,rr,zz,phi,Ht,zot,bt)
    theta =3*pi/2
    rr = 1 + rpt*np.cos(theta)
    zz = phi/omet + rpt*np.sin(theta)
    y4 = G_function_spiral(Fo,rr,zz,phi,Ht,zot,bt)
    th2 = (y1+y2+y3+y4)/4
    return th2


def G_function_fls(Fo, Ht = 1/0.0005,zot = 80):
    if Fo == 0:
        G = 0.0
    else:
        rbb = 1/Ht
        zob = zot*rbb
        tt = 9*rbb*rbb*Fo
        G = fct_fls(tt,rbb,zob)/(2*pi)
    return G
def G_function_fls_loc(Fo,r,z,h,zo):
    def g_integrandz(x,Fo,r,z,h,zo):
        rp = np.sqrt(r**2+(z-x)**2)
        rm = np.sqrt(r**2+(z+x)**2)
        y = erfc(rp/(2*np.sqrt(Fo)))/rp-erfc(rm/(2*np.sqrt(Fo)))/rm
        return y
    if Fo == 0:
        G = 0.0
    else:
        I = quad(g_integrandz, zo, zo+h, args=(Fo,r,z,h,zo))
        G = I[0]/(4*pi)
    return G

def  G_function(Fo):
#
# fonction cylindrique (Cooper 1976)
#
    Ca=1.128379;C0=-0.5;C1=0.2756227;C2=-0.1499385;C3=0.0617932;C4=-0.01508767;C5=0.001566857
    euler = 0.57721566490153286060651209008240243104215933593992
    if Fo <= 6.124633:
        sFo = np.sqrt(Fo)
        G = sFo/(2*pi)*(Ca+C0*sFo+C1*Fo+C2*Fo**(1.5)+C3*Fo**(2.0)+C4*Fo**(2.5)+C5*Fo**(3.0))
    else:
        z = np.log(4*Fo/np.exp(euler))
        G = (2*z*(8*Fo*(1+2*Fo)-1.0-3*z)+16*Fo+pi**2+3)/(128.0*pi*Fo**2)
    return G

def  G_function_bernier(Fo , p=1):
#
# fonction cylindrique (INgersoll 1954)
#
    A = np.array([[-0.89129,0.36081,-0.05508,3.59617e-3,0,0], \
        [-1.454099,0.8993336,-0.311928,0.061119,-0.00478046,0], \
        [-3.007691,2.256059,-0.7928093,0.134293,-0.00858244,0], \
        [-9.141771,11.7025,-7.09574,2.269837,-0.3669166,0.023587]])
    if p == 1:
        i = 0
    elif p ==2:
        i = 1
    elif p == 5:
        i = 2
    elif p == 10:
        i = 3
    else:
        print ('le deuxieme argument doit être 1,2,5 ou 10')
        G = -999
        return G
    c = A[i,:]
    x = np.log10(Fo)
    arg = c[0]+c[1]*x+c[2]*x**2+c[3]*x**3+c[4]*x**4+c[5]*x**5
    G = 10**arg
    return G


def g_function_fls(tt,rbb = 0.0005,zob = 0.04):
    if tt == 0:
        g = 0.0
    else:
        g = fct_fls_claesson(tt,rbb,zob)
    return g

def g_function_claessonb(Fo,Ht=2000.0):
    if Fo == 0:
        G = 0.0
    else:
        G = fct_fls_claessonb(Fo,Ht)/(4*pi*Ht)
    return G

def h_function_claesson(Fo,Ht=2000.0):
    if Fo == 0:
        G = 0.0
    else:
        G = (Fo*fct_fls_claessonb(Fo,Ht) - fct_fls_claesson_n(Fo,Ht))/(4*pi*Ht)
    return G

def rinfc(Fo,xi):
    def fct(x):
        y = Fo - (0.0649*x**2+0.2409-0.5*np.log(x)-(1.0/(2.0*x**2)-1.0/(3.0*x**4)+2.0/(9.0*x**6)-2.0/(15.0*x**8)))
        return y
    try:
        z = newton(fct,xi)
    except RuntimeError:
        print ("error")
        z = 0.1
    return z

def G_Hart(rinft):

    if rinft>3:
        G = 0.5*(np.log(rinft)-0.9818+2.0/rinft**2-2.0/rinft**4)/pi
    else:
        x = 4.0/rinft**2
        G = (np.log(rinft)-0.9818+2.0/rinft**2-2.0/rinft**4)
        s = -1
        i = 3
        ok = False
        while (ok==0):
            s = -s
            delG =  s*x**i/(2*i*factorial(i))
            G = G+delG
            if(abs(delG) < 1e-6):
                ok = True
            else:
                i = i+1
        G = G/(2*pi)
    return G

def  g_function_hart_couv(Fo,rt=1):
#
# Fonction grand - G évaluée approximée par la ligne source de Kelvin de
# Hart Couvillon appliquée à une conduite (Chapitre 3 du livre H_C)
# Fo = al*t/rb^2;
# rt = r/rb
#
    if Fo < 0.2:
        G = 0
        rinf = 0
        return G,rinf
    rinfi = 4*np.sqrt(Fo)/rt    #ceci est (rinf/r)
    #G = pwrfor(1/rinfi)/(2*pi);
    if rinfi < 15:
        rinf = rinfc(Fo,rinfi)
    else:
        rinf = rinfi
    G = G_Hart(rinf)
    return G,rinf


def G_function_ils(Fo,rt=1):
    if Fo == 0:
        G = 0.0
    else:
        z = rt*rt/(4.0*Fo)
        G = expn(1,z)/(4.0*pi)
    return G

def g_function_ilsn(Fo,rt=1):
    if Fo == 0:
        G = 0.0
    else:
        z = rt*rt/(4.0*Fo)
        G = Fo*(expn(1,z)-expn(2,z))/(4.0*pi)
    return G


def fundb(x,ht):
    y = np.exp(-x*x)*(4.0*ierf(x*ht)-ierf(2.0*x*ht))/(x*x)
    return y

def fun_n(x,ht):
    y = np.exp(-x*x)*(4.0*ierf(x*ht)-ierf(2.0*x*ht))/(4*x*x*x*x)
    return y

def fct_fls_claesson(t,r=0.0005,a=0.04):
    def fund(x,rt,dt):
        y = np.exp(-x*x*rt*rt)*(2.0*ierf(x)+2.0*ierf((1.0+2.0*dt)*x)-ierf(2.0*x*(1.0+dt))-ierf(2.0*dt*x))/(x*x)
        return y
    xi = 3.0/(2.0*np.sqrt(t))
    y2 = quad(fund,xi,Inf,args=(r,a))[0]/2.0
    return y2

def fct_fls_claessonb(t,ht):
    xi = 1.0/(2.0*np.sqrt(t))
    y2 = quad(fundb,xi,Inf,args=(ht))[0]
    return y2

def fct_fls_claesson_n(t,ht):
    xi = 1.0/(2.0*np.sqrt(t))
    y2 = quad(fun_n,xi,Inf,args=(ht))[0]
    return y2

def I_function(xx):
    z = xx*xx
    I = expn(1,z)/(2.0)
    return I


def W_function(u):
    I = expn(1,u)
    return I

def sdr_pipe(pipe,sdr):

# facteurs de conversion
# D1 est le diametre int?rieur en m?tres
# D2 est le diametre ext?rieur en m?tres

    fac5 = 0.3048            # m/pied

    if pipe == 0.5:
        D2 = 0.84/12*fac5
    elif pipe == 0.75:
        D2 = 1.05/12*fac5
    elif pipe == 1:
        D2 = 1.315/12*fac5
    elif pipe == 1.25:
        D2 = 1.66/12*fac5
    elif pipe == 1.5:
        D2 = 1.9/12*fac5
    elif pipe == 2:
        D2 = 2.375/12*fac5
    elif pipe == 3:
        D2 = 3.5/12*fac5
    elif pipe == 4:
        D2 = 4.5/12*fac5
    elif pipe == 5:
        D2 = 5.562/12*fac5
    elif pipe == 6:
        D2 = 6.625/12*fac5
    else:
        print ('pipe number not defined'  )
        D2 = 0
    t = D2/sdr
    D1 = D2 - 2*t
    return D1,D2


def  Rb_linesource(kg,ks,rb,rp,xc):

#
# kg : grout
# ks : sol
# rb : rayon puits
# rp : rayon ext?rieur
# xc : demi-distance des tuyaux
    sig = (kg-ks)/(kg+ks)
    b11 = xc/rb
    b22 = b11
    b12 = 2*xc/rb
    la1 = rb/rp
    la2 = rb/xc
    la3 = la2/(2*la1)
    R11 = (np.log(la1)+np.log(la2/2)+sig*np.log(la2**4/(la2**4-1)))/(2*pi*kg)
    Ra = (np.log(2*xc/rp)+sig*np.log((1+b11**2)/(1-b11**2)))/(pi*kg)
    Rb = R11/2
    return Rb,Ra


def  Rb_linesource_2UT(kg,ks,rb,rp,xc,cas = '13-24'):

#
# kg : grout
# ks : sol
# rb : rayon puits
# r1 : rayon ext?rieur
# xc : demi-distance des tuyaux
# Rp : R?sistance des tuyaux
    sig = (kg-ks)/(kg+ks)
    xcp = xc/rb
    r1p = rp/rb
    Rb = -(np.log(4*r1p*xcp**3) + sig*np.log(1 - xcp**8))/(8*pi*kg)
    if cas == '13-24':
        Ra = (np.log(2*xcp/r1p) + sig*np.log((1+xcp**2)/(1-xcp**2)))/(2*pi*kg)
    else:
        Ra = (np.log(xcp/r1p) + sig*np.log((1+xcp**4)/(1-xcp**4)))/(2*pi*kg)
    return Rb,Ra




def  Rb_Paul(kg,rb,rp,cas):
#
# calcul resistance Redmund
# cas = upper(cas)
    if cas == 'a' or cas == 'A':
        beta0 = 20.10
        beta1 = -0.94447
    elif cas == 'b' or cas == 'B':
        beta0 = 17.44
        beta1 = -0.6052
    elif cas == 'c' or cas == 'C':
        beta0 = 21.90587
        beta1 = -0.3796
    else:
        print ('error, cas = a b or c' )
    Sb = beta0*(rb/rp)**beta1
    Rg = 1/(kg*Sb)
    return Rg


def  Rb_Sharqawi(kg,rb,rp,xc):
#
# calcul resistance Redmund
# cas = upper(cas)
    Rb = 1/(2*pi*kg)*(-1.49*(xc/rb)+0.656*np.log(rb/rp)+0.436)
    return Rb

def  calcul_Ra_shape_factor(kg,rp,xc):
#
    Sb = 2*pi/np.arccosh(2*(xc/rp)**2-1)
    Ra = 1/(kg*Sb)
    return Ra




def  calcul_mois_jour(jours):
# fonction donnant le jour et le mois de l'annee correspondant au nombre de jours depuis le d?but de l'ann?e
# [mois,jour] = calcul_jour_mois(100)
# mois = 4, jour = 10
    jrr = array([31,59,90,120,151,181,212,243,273,304,334,365])
    mois = cherche_index(jours,jrr)   # im donne le mois de l'annee
    if (mois>0):
        jour = jours - jrr[mois]
    else:
        jour = jours
    return mois+1,int(jour)

def  calcul_mois_jour_heure(heures):
# fonction donnant le jour et le mois de l'annee correspondant au nombre d'herues depuis le d?but de l'ann?e
# mois,jour,heure = calcul_mois_jour_heure(100)
# mois = 4, jour = 10
    hrr = array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])
    mois = cherche_index(heures,hrr)   # im donne le mois de l'annee
    if (mois>0):
        reste = heures - hrr[mois-1]
    else:
        reste  = heures
    jour = ceil((reste+1)/24.0)
    heure = reste -(jour-1)*24
    return mois+1,int(jour),int(heure)


def  Rp_fct(mp,ro,ri,kp,Trefk,fluid):
    # Fonction calculant la r?sistance de convection et de conduction dans le
    # tube de HDPE dans le cas ou le fluide caloporteur est de l'eau
    # mp : debit en kg/sc
    # Tref en celcius
    # ro : rayon exterieur du tuyau
    # rayon interieur du tuyau
    # kp : conductivite du tuyau
    patm = 101.325*1000.0
    di = 2*ri
    visc = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
    Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
    kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
    rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
    Re = (4*mp)/(pi*di*visc)           # nombre de Reynold
    cas = 1
    if (Re>2300):
        if cas == 1:
                # Gnielienski
                f = (0.79*np.log(Re)- 1.64)**-2
                Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
                # Dittus-Boetter
        else:
                Nu1 = 0.023*Re**(0.8)*Pr**0.4
                Nu2 = 0.023*Re**(0.8)*Pr**0.3
                Nud = (Nu1+Nu2)/2
    else:
        Nud = 3.6
        print('laminar')
    # calcul du coefficient de convection
    hf = (Nud*kf)/(di)
    # r?sistance de convection
    rconv = 1/(hf*pi*di)
    # r?sistance de conduction
    rcond = np.log(ro/ri)/(2.0*pi*kp)
    Rp = rcond + rconv
    return Rp,rcond,rconv



def Rb_multipole(kb,k,rb,rp,Rp,Jp,z):

    #
    # kb : grout
    # k : sol
    from copy import  deepcopy
    beta = 2*pi*kb*Rp
    sig = (kb-k)/(kb+k)
    la1 = rb/rp
    N = len(z)
    Ro = np.zeros((N,N))
    r = abs(z)
    Nm = N*(N+1)/2.0
    for i in range (0,N):
        Ro[i,i] = (beta+np.log(la1)+sig*np.log(rb**2/(rb**2-r[i]**2)))/(2*pi*kb)
    for i in range(0,N-1):
        for j in range (i+1,N):
            rmn = abs(z[i]-z[j])
            Ro[i,j] = (np.log(rb/rmn)+sig*np.log(rb**2/np.abs(rb**2-np.conj(z[i])*z[j])))/(2*pi*kb)
            Ro[j,i] = Ro[i,j]
    if Jp == 0:
        R = Ro
    else:
        P = np.zeros((N,Jp)).astype(complex)
        Pn = np.zeros((N,Jp)).astype(complex)
        F = np.zeros((N,Jp)).astype(complex)
        R = np.zeros((N,N)).astype(complex)
        for m in range (0,N):
            q = np.zeros(N)
            q[m] = 1
            ok = 0
            compt = 1
            compt_max = 100
            delta = 1e-5
            while (ok==0):
                for i in range (0,N):
                    for ik in range (0,Jp):
                        zz = Calcul_F(q,P,rp,rb,kb,sig,i,ik,Jp,N,z)
                        F[i,ik] = zz
                        Pn[i,ik] = (-1+(ik+1)*beta)/(1+(ik+1)*beta)*np.conj(zz)
                d = Pn-P
                dell = np.abs(d)
                ep = dell.max()
                if ep < delta:
                    ok = 1
                else:
                    P = deepcopy(Pn)
                    compt = compt + 1
                    if compt > compt_max:
                        err = 1
                        ok = 1
            for im in range (m,N):
                s1 = 0
                s2 = 0
                for i in range(0,N):
                    if i != im:
                        for ij in range (0,Jp):
                            s1 = s1 + P[i,ij]*rp**(ij+1)/(z[im] -z[i])**(ij+1)
                    for ij in range (0,Jp):
                         s2 = s2 + P[i,ij]*rp**(ij+1)*np.conj(z[im])**(ij+1)/(rb**2-np.conj(z[im])*z[i])**(ij+1)
                R[m,im] = Ro[m,im] + np.real(s1+sig*s2)
                R[im,m] = R[m,im]
    Rm = np.linalg.inv(R)
    K = np.zeros((N,N))
    for i in range (0,N):
        K[i,i] = sum(np.real(Rm[i,:]))
    for i in range (0,N-1):
        for j in range (i+1,N):
            K[i,j] = -np.real(Rm[i,j])
            K[j,i] = K[i,j]
    Kb = sum(np.diag(K))
    Rb = 1/Kb
    if N == 2:
        Ka = (K[0,0] + 2*K[0,1])/2.0
        Ra = 1/Ka
    elif N==4:
        Ka1 = (K[0,0] + 2*K[0,1]+2*K[0,2])
        Ka2 = (K[0,0] + 2*K[0,1]+2*K[0,3])
        Ra = 1/Ka1
    return Rb,Ra


def  Calcul_F(q,P,rp,rb,kb,sig,m,km,Jp,N,z):

    k = km+1
    sa = 0
    for i in range (0,N):
        if i != m:
            den = 2*pi*kb*k*(z[i]-z[m])**k
            sa = sa + q[i]*rp**k/den
    sc = 0
    for i in range (0,N):
        if i != m:
            for ij in range(0,Jp):
                num = P[i,ij]*comb(ij+k,ij)*rp**(ij+1)*(-rp)**k
                den = (z[m] - z[i])**(ij+k+1)
                sc = sc+num/den
    sb = 0
    for i in range (0,N):
        num = q[i]*rp**k*np.conj(z[i])**k
        den = 2*pi*kb*k*(rb**2-z[m]*np.conj(z[i]))**k
        sb  = sb + num/den
    sd = 0
    for i in range (0,N):
        for ij in range (0,Jp):
            nj = min(ij+1,k)
            sj = 0
            for jp in range (0,nj+1):
                num1 = comb(ij+1,jp)*comb(ij+k-jp,ij)
                num2 = rp**(ij+k+1)*z[m]**(ij+1-jp)*np.conj(z[i])**(k-jp)
                den = (rb**2-z[m]*np.conj(z[i]))**(k+ij+1-jp)
                sj = sj+num1*num2/den
            sd = sd +np.conj(P[i,ij])*sj
    F = sa + sig*(sb+sd) + sc
    return F

def  Tp_ashrae(nx,ny,d,qp,al,ks,n_annees,nrings=3,H = -999):

    # nx : nombre de puits en x (colonnes)
    # ny : nombre de puits en y (rang?es)
    # d : distance entre les puits (On suppose toujours la m?me )_
    # al : diffusivité du sol (m2/jour)
    # ks : conductivité du sol (W/mK)
    # n_annees : nombre d'ann?es
    # na : nombre d'anneau
    if nx*ny ==1:
        return 0
    ri = d/2.0
    dr = 1.5
    tf = n_annees*365     # temps final en jours
    r1 = ri
    Tp = 0
    r1 = ri
    for i in range (0,nrings):
        r2 = r1+dr
        rm = (r1+r2)/2.0
        X  = rm/(2.0*np.sqrt(al*tf))
        Ix  = I_function(X)
        Tm = qp/(2.0*pi*ks)*Ix
        Tp = Tp+Tm*pi*(r2**2-r1**2)/d**2
        r1 = r2
    # Formule exacte pour Tp non corrig?e
    z = ri**2/(4.0*al*tf)
    E2 = np.exp(-z) - z*expn(1,z)
    q2 = qp*tf*E2
    CC = ks/al
    Tpf = q2/d**2/CC
    # correction pour les effets de bout
    ntot = nx*ny # nombre d puits
    if nx==1 :
        n1 = min(2,ntot)    # nombre coins
        n2 = max(ny-n1,0)    # nombre de puits sur la p?riph?rie
        n3 = 0
        n4 = 0
    elif ny==1:
        n1 = min(2,ntot)    # nombre coins
        n2 = max(nx-n1,0)    # nombre de puits sur la p?riph?rie
        n4 = 0
        n3 = 0
    else:
        n1 = 0
        n2 = 4
        n3 = 2*(nx-2)+2*(ny-2)    # nombre de puits sur la p?riph?rie
        n4 = ntot - n3 - n2        # nombre de puits int?rieurs
    #Tpn = (n4+0.5*n3+0.25*n2+0.1*n1)/ntot*Tp
    Wf = (nx-1)*d
    Lf = (ny-1)*d
    if H > 0:
        Cf = (H*2*(Wf+Lf)+Wf*Lf)/(H*2*(Wf+Lf))
    else:
        Cf = 1
    Tpn = (n4+0.75*n3+0.5*n2+0.25*n1)/(ntot*Cf)*Tp
    return Tpn


def Tp_ils(nx,ny,d,qp,al,ks,n_annees):

    # nx : nombre de puits en x (colonnes)
    # ny : nombre de puits en y (rang?es)
    # d : distance entre les puits (On suppose toujours la m?me )_
    # al : diffusivit? du sol (m2/jour)
    # ks : conductivit? du sol (W/mK)
    # n_annees : nombre d'ann?es

    ntot = nx*ny
    Tpp = np.zeros((ntot,1))
    x = np.zeros((ntot,1))
    y = np.zeros((ntot,1))
    tf = n_annees*365     # temps final en jours
    k = 0
    dx = d
    dy = d
    for i in range(0,nx):
        for j in range(0,ny):
            x[k] = dx*(i-1)
            y[k] = dy*(j-1)
            k = k+1
    s1 = 0
    for i in range (0,ntot-1):
        for j in range(i+1,ntot):
            dist = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
            xx  = dist[0]/(2*np.sqrt(al*tf))
            Ixx  = I_function(xx)
            s1 = s1 + qp/(2*pi*ks)*Ixx
    Tp_nouv = 2*s1/ntot
    return Tp_nouv


def Tp_flsn(nx,ny,d,qp,al,ks,n_annees,H=100.0,zob = 0.04):

    # nx : nombre de puits en x (colonnes)
    # ny : nombre de puits en y (rang?es)
    # d : distance entre les puits (On suppose toujours la m?me )_
    # al : diffusivit? du sol (m2/jour)
    # ks : conductivit? du sol (W/mK)
    # n_annees : nombre d'ann?es
    ts = H**2/(9*al)         #temps caract?ristique Eskilson en jour
    dx = d/H
    dy = d/H
    ntot = nx*ny
    Tpp = np.zeros((ntot,1))
    x = np.zeros((ntot,1))
    y = np.zeros((ntot,1))
    gbb = np.zeros((ntot,ntot))
    tf = n_annees*365/ts     # temps final en jours
    zt = np.zeros(ntot).astype(complex)
    k = 0
    for i1  in range(0,nx):
        for i2 in range(0,ny):
            zt[k] = i1*dx+1j*dy*i2
            k = k+1
    s1 = 0
    for i in range (0,ntot-1):
        for j in range(i+1,ntot):
            rb = np.abs(zt[j]-zt[i])
            Ixx  = fct_fls(tf,rb,zob)
            s1 = s1 + Ixx
    Tp_nouv = s1*qp/(pi*ks*ntot)
    return Tp_nouv

def Tp_eight(nx,ny,d,qp,al,ks,n_years,H=100.0,conf = 'R'):

    # d : distance entre les puits (On suppose toujours la m?me )_
    # alhr : diffusivit? du sol (m2/jour)
    # ks : conductivit? du sol (W/mK)
    #
    tf = n_years*8760
    alhr = al/24
    Href = 100
    bh = d/H
    fact = Href/H
    ntot = nx*ny
    if nx==1 :
        n1 = min(2,ntot)    # nombre coins
        n2 = max(ny-n1,0)    # nombre de puits sur la p?riph?rie
        n3 = 0
        n4 = 0
    elif ny==1:
        n1 = min(2,ntot)    # nombre coins
        n2 = max(nx-n1,0)    # nombre de puits sur la p?riph?rie
        n4 = 0
        n3 = 0
    else:
        n1 = 0
        n2 = 4
        n3 = 2*(nx-2)+2*(ny-2)    # nombre de puits sur la p?riph?rie
        n4 = ntot - n3 - n2        # nombre de puits int?rieurs
    Fo1 = alhr*tf/d**2*fact
    Fo2 = alhr*tf/(2*d**2)*fact
    x1 = 1/(4*Fo1)
    x2 = 1/(4*Fo2)
    xx1 = np.sqrt(x1)
    xx2 = np.sqrt(x2)
    Th8 = qp/(pi*ks)*(expn(1,x1)+expn(1,x2))
    if conf == 'R':
        a = 1.95005 + 0.105215/bh - 55.6543*bh**2
        b = 0.28
        c = 0.45
        d = 0
    else:
        a = 0
        b = 0.95
        c = -0.28174*np.log(bh) - 0.23546
        d = 0.05
    Tp8 = Th8*(a*n4 + b*n3 + c*n2 + d*n1)/ntot
    return Tp8

def Tp_Capozza(nx,ny,d,qp,al,ks,n_years,H=100.0):

    # d : distance entre les puits (On suppose toujours la m?me )_
    # alhr : diffusivit? du sol (m2/jour)
    # ks : conductivit? du sol (W/mK)
    # tf
    bh = d/H
    tss = H**2/(9*al)
    tf = n_years*365
    tt = tf/tss
    beta = 6/np.log(0.07)
    x = np.log(tt)
    nb = nx*ny
    Nbv = np.array([1,4,8,16,25,100,225])
    bv = np.array([0,1.03,1.26,1.52,1.59,2.66,2.32])
    fv = np.array([1,1.24,1.245,1.25,1.30,1.31,1.41])
    b = np.interp(nb,Nbv,bv)
    f = np.interp(nb,Nbv,fv)
    ratio = 1 +(b + (1+b)*(f-1)/0.05*(bh-0.05))*np.exp(-(x-3)/beta)
    Tpi = Tp_ils(nx,ny,d,qp,al,ks,n_years)
    Tp = Tpi/ratio
    return Tp

def Tp_fls(nx,ny,d,qp,al,ks,n_annees,H=100.0,zob = 0.04):

    # nx : nombre de puits en x (colonnes)
    # ny : nombre de puits en y (rang?es)
    # d : distance entre les puits (On suppose toujours la m?me )_
    # al : diffusivit? du sol (m2/jour)
    # ks : conductivit? du sol (W/mK)
    # n_annees : nombre d'ann?es
    ts = H**2/(9*al)         #temps caract?ristique Eskilson en jour
    dx = d/H
    dy = d/H
    ntot = nx*ny
    Tpp = np.zeros((ntot,1))
    x = np.zeros((ntot,1))
    y = np.zeros((ntot,1))
    gbb = np.zeros((ntot,ntot))
    tf = n_annees*365/ts     # temps final en jours
    zt = np.zeros(ntot).astype(complex)
    k = 0
    for i1  in range(0,nx):
        for i2 in range(0,ny):
            zt[k] = i1*dx+1j*dy*i2
            k = k+1

    for i in range (0,ntot):
#        gbb[i,i] = 0
        for j in range(i+1,ntot):
            if(i!=j):
                rb = np.abs(zt[j]-zt[i])
                Ixx  = fct_fls(tf,rb,zob)
                gbb[i,j] = Ixx
                gbb[j,i] = Ixx
        Tpp[i] = np.sum(gbb[i,:])
    Tp_nouv = np.mean(Tpp)*qp/(2.0*pi*ks)
    return Tp_nouv



def Tp_Bernier(nx,ny,d,qp,al,ks,n_annees,H =100.0):
    b = np.array([7.8189,-6.4270E+01,+1.5387E+02,-8.4809E+01,+3.4610E+00,-9.4753E-01, \
        -6.0416E-02,+1.5631E+00,-8.9416E-03,+1.9061E-05,-2.2890E+00, \
        +1.0187E-01,+6.5690E-03,-4.0918E+01,+1.5557E+01,-1.9107E+01, \
        +1.0529E-01,+2.5501E+01,-2.1177E+00,+7.7529E+01,-5.0454E+01, \
        +7.6352E+01,-5.3719E-01,-1.3200E+02,+1.2878E+01,+1.2697E-01, \
        -4.0284E-04,-7.2065E-02,+9.5184E-04,-2.4167E-02,+9.6811E-05, \
        +2.8317E-02,-1.0905E-03,+1.2207E-01,-7.1050E-03,-1.1129E-03,-4.5566E-04])
    x = nx*d
    y = ny*d
    zp = max(x,y)
    zm = min(x,y)
    A = zp/zm
    ts = H**2/(9*al)         #temps caract?ristique Eskilson en jour
    dx = d/H
    tf = n_annees*365/ts     # temps final en jours
    xf = np.log(tf)
    nb = nx*ny
    sum = 0
    for i in range(0,4):
        sum = sum + b[i]*dx**i
    for i in range(0,3):
        j = i+4
        sum = sum + b[j]*xf**(i+1)
    for i in range(0,3):
        j = i+7
        sum = sum + b[j]*nb**(i+1)
    for i in range(0,3):
        j = i+10
        sum = sum + b[j]*A**(i+1)
    for i in range(0,2):
        j = i + 13
        sum = sum + b[j]*dx*xf**(i+1)
    for i in range(0,2):
        j = i + 15
        sum = sum + b[j]*dx*nb**(i+1)
    for i in range(0,2):
        j = i + 17
        sum = sum + b[j]*dx*A**(i+1)
    for i in range(0,2):
        j = i + 19
        sum = sum + b[j]*dx*dx*xf**(i+1)
    for i in range(0,2):
        j = i + 21
        sum = sum + b[j]*dx*dx*nb**(i+1)
    for i in range(0,2):
        j = i + 23
        sum = sum + b[j]*dx*dx*A**(i+1)
    for i in range(0,2):
        j = i + 25
        sum = sum + b[j]*xf*nb**(i+1)
    for i in range(0,2):
        j = i + 27
        sum = sum + b[j]*xf*A**(i+1)
    for i in range(0,2):
        j = i + 29
        sum = sum + b[j]*xf*xf*nb**(i+1)
    for i in range(0,2):
        j = i + 31
        sum = sum + b[j]*xf*xf*A**(i+1)
    for i in range(0,2):
        j = i + 33
        sum = sum + b[j]*nb*A**(i+1)
    for i in range(0,2):
        j = i + 35
        sum = sum + b[j]*nb*nb*A**(i+1)
    Tp_nouv = sum*qp/(2.0*pi*ks)
    return Tp_nouv

def  funa(x,r,beta):
        return erfc(beta*x)/np.sqrt(x**2-r**2)

def fct_hij(t,rr,z1,z2,h1,h2):

    xi = 3/(2*np.sqrt(t))

    def fct_massimo(x):
        x1 = (z2 - z1 + h2)*x
        x2 = (z2 - z1)*x
        x3 = (z2 - z1 - h1)*x
        x4 = (z2 - z1 + h2 - h1)*x
        x5 = (z2 + z1 + h2)*x
        x6 = (z2 + z1)*x
        x7 = (z2 + z1 + h1)*x
        x8 = (z2 + z1 + h2 + h1)*x
        y = np.exp(-x*x*rr*rr)*(ierf(x1)- ierf(x2) + ierf(x3) - ierf(x4)+ierf(x5)- ierf(x6) + ierf(x7) - ierf(x8))/(x*x)
        return y

    if t == 0:
        h = 0
    else:
        h = quad(fct_massimo,xi,Inf)[0]
    return h

def fct_fls(t,r=0.0005,zob=0.04):
#
# fonction g_minuscule calculé par la ligne source allant de d ? H
# r  = rb/H
# t = t/ts
# a = d/H
    a = zob
    epp = np.spacing(1)
    if t == 0:
        y = 0
    else:
        beta = 3.0/(2.0*np.sqrt(t))
        am = np.sqrt(r**2+1);
        a1 = np.sqrt(r**2+(1+2*a)**2)
        a2 = np.sqrt(r**2+(2+2*a)**2)
        a3 = np.sqrt(r**2+4*a**2)
        D1 = am*erfc(beta*am)-r*erfc(beta*r)-(np.exp(-beta**2*am**2)-np.exp(-beta**2*r**2))/(beta*np.sqrt(pi))
        y1 = quad(funa,r+epp,am,args=(r,beta))[0]
        A = y1-D1
        D2 = a1*erfc(beta*a1)-(a3*erfc(beta*a3)+a2*erfc(beta*a2))/2.0-(np.exp(-beta**2*a1**2)-(np.exp(-beta**2*a3**2)+np.exp(-beta**2*a2**2))/2.0)/(beta*np.sqrt(pi))
        y2b = quad(funa,a1,a2,args=(r,beta))[0]
        y3b = quad(funa,a3,a1,args=(r,beta))[0]
        B = D2+y2b*(1+a)-a*y3b
        y = A-B
    return y




def compute_g_function(zb,t,rbb=0.0005,zob =0.04):
#
# calul de la valeur  de la fontion g COMMENCANT a UNE PROFIONDEUR D
# zb : liste des puits , coordonn?es sous forme de nombre  complexe z =
# x+yj
# rr = rb/H
# t = t/ts
    nb = len(zb)
    gb = np.zeros(nb)
    gbb = np.zeros((nb,nb))
    g = fct_fls(t,rbb,zob)
    for jj in range(0,nb):
        gbb[jj,jj] = g
        for k in range(jj+1,nb):
            if(k!=jj):
                rb = np.linalg.norm(zb[k]-zb[jj])
                val = fct_fls(t,rb,zob)
                gbb[jj,k]= val
                gbb[k,jj] = val
        gb[jj] = sum(gbb[jj,:])
    g2 = np.mean(gb)
    return g2


def Pulses_ashrae(q_sol,nbloc = 4,flag = 1):
#
# calcul les pulse qa, qm_ch,qm_cl,qh_ch,q_cl, à partir des charges horaires au sol
# nbloc : bloc horaire , par défaut = 4
# flag , 1 , mois réels (défaut), 2 mois = 730 heures
#
    qa = np.mean(q_sol)
    qmois = np.zeros(12)
    #
    # calcul des charges mensuelles (puissance moyenne  et puissance maximale)
    #
    i1 = 0
    for i in range (0,12):      # sommation pour chacun des mois
        if flag == 1:
            i2 = i1+hrm[i]                    # i1 première heure du mois, i2, dernière heure du mois
        else:
            i2 = i1+ 730
        qmois[i] = np.mean(q_sol[i1:i2])    # puissance moyenne chauffage
        i1 = i2
    qm_ch = max(qmois)       # Pulses mensuels chaufffage pour chaque zone
    qm_cl = min(qmois)       # Pulses mensuels climatisation pour chaque zone
    qm_ch = max(qm_ch,0)
    qm_cl = min(qm_cl,0)
    #
    # Calcul des moyennes sur bloc horaire
    #
    ibloc = int(8760/nbloc)   # nombre de blocs par années
    q_bloc = np.zeros(ibloc)       # pulses bloc horaire
    i1 = 0
    for i  in range(0,ibloc):
        i2 = int((i+1)*nbloc)
        q_bloc[i] = np.mean(q_sol[i1:i2])    # puissance moyenne chauffage
        i1 = i2
    qh_ch = max(q_bloc)       # Pulses mensuels chaufffage pour chaque zone
    qh_cl = min(q_bloc)       # Pulses mensuels climatisation pour chaque zone
    qh_ch = max(qh_ch,0)
    qh_cl = min(qh_cl,0)
    return qa,qm_ch,qm_cl,qh_ch,qh_cl

def Delta_T_earth(x,ts,tshift,al,amp):

    # al  en m2/jour
    #  x en  m
    ome = 2.0*pi/365.0
    beta = np.sqrt(ome/(2.0*al))
    a = np.sqrt(pi/(365.0*al))
    b = ts - tshift - x/2.0*np.sqrt(365.0/(pi*al))
    dT =  -amp*np.exp(-x*beta)*np.cos(ome*b)
    dTmax =  amp*np.exp(-x*beta)
    return dT,dTmax

def Calcul_R_horizontal_1puits(X,zt,ks):

    Xi = 2.0*zt*X
    R = (I_function(X)-I_function(Xi))/(2*pi*ks)
    return R


def Calcul_Rg_horizontal_1puits(X,zt,ks):

    Fo = 1.0/(4*X*X)
    xt = 2.0*zt
    R = (G_function(Fo)-G_function_int(Fo,xt))/(ks)
    return R

def Calcul_R_horizontal_2puits(X1,zt1,zt2,ks):


    X2 = (zt2-zt1)*X1
    X3 = 2*zt1*X1
    X4 = (zt1+zt2)*X1
    R1 = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    X3 = 2*zt2*X1
    R2 = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    R = (R1+R2)/2.0
    return R

def Calcul_R_horizontal_4puits(X1,zt1,zt2,dt2,ks):

    d1 = zt2 - zt1
    X2 = d1*X1
    X3 = dt2*X1
    X4 = sqrt(dt2**2+d1**2)*X1
    X5 = 2*zt2*X1
    X6 = (zt2+zt1)*X1
    X7 = sqrt(dt2**2+(2*zt2)**2)*X1
    X8 = sqrt(dt2**2+(zt1+zt2)**2)*X1
    Rs1 = ((I_function(X1)+I_function(X2)+I_function(X3)+I_function(X4))-(I_function(X5)+I_function(X6)+I_function(X7)+I_function(X8)))/(2*pi*ks)     # mK/W
    X5 = 2*zt1*X1
    X7 = sqrt(dt2**2+(2*zt1)**2)*X1
    Rs2 = ((I_function(X1)+I_function(X2)+I_function(X3)+I_function(X4))-(I_function(X5)+I_function(X6)+I_function(X7)+I_function(X8)))/(2*pi*ks)     # mK/W
    R = (Rs1+Rs2)/2.0
    return R

def  Compute_Taxial(Ra,Rb,CCC,L):
    eta = L/(CCC*np.sqrt(Ra*Rb))
    xsi = np.sqrt(Ra/(4.0*Rb))
    dz = L/15.0
    z = np.arange(0,L+dz,dz)
    n = len(z)
    Td = np.zeros(n)
    Tu = np.zeros(n)
    nn = 2*n
    Ty = np.zeros(nn)
    zy = np.zeros(nn)

    xin = np.sqrt(Ra/Rb)/2.0
    for i in range(0,n):
        zx = z[i]/L
        x = eta*z[i]/L
        y = eta*(1 - zx)
        den = np.cosh(eta) + xsi*np.sinh(eta);
        Td[i] = (np.cosh(y) + xsi*np.sinh(y))/den
        Tu[i] = (np.cosh(y) - xsi*np.sinh(y))/den
        Ty[i] = Td[i]
        Ty[nn-i - 1] = Tu[i]
        zy[i] = z[i]
        zy[nn-i-1] = z[i]
    return z,Td,Tu,zy,Ty


def  Calcul_TBeir(R12,R1,CCC,L,Rs):
    N12 = L/(CCC*R12)
    Ns1 = L/(CCC*(R1 + Rs))
    Ns2 = Ns1
    num = np.sqrt((Ns1 - Ns2)**2 + 4.0*((N12+Ns1)*(N12+Ns2) - N12**2))
    a1 = (-(Ns1-Ns2) + num)/2.0
    a2 = (-(Ns1-Ns2) - num)/2.0
    C1 = (Ns1 + a2)*np.exp(a2)/((Ns1+a2)*np.exp(a2) - (Ns1 + a1)*np.exp(a1))
    C2 = 1 - C1
    C3 = (Ns1 + N12 + a1)*C1/N12
    C4 = (Ns1 + N12 + a2)*C2/N12
    C5 = C1*(1 + Ns2/Ns1*(N12 + Ns1 + a1)/N12)*(np.exp(a1)-1)/a1
    C6 = (1-C1)*(1 + Ns2/Ns1*(N12 + Ns1 + a2)/N12)*(np.exp(a2)-1)/a2
    Q = C5 + C6
    dz = L/15.0
    z = np.arange(0,L+dz,dz)
    n = np.len(z)
    Td = np.zeros(n)
    Tu = np.zeros(n)
    Tu2 = np.zeros(n)
    nn = 2*n
    Ty = np.zeros(nn)
    Ty2 = np.zeros(nn)
    zy = np.zeros(nn)
    for i in range(0,n):
        x1 = a1*z[i]/L;
        x2 = a2*z[i]/L;
        x3 = a1*(2*L - z[i])/L
        x4 = a2*(2*L - z[i])/L
        Td[i] = C1*np.exp(x1) + C2*np.exp(x2);
        Tu[i] = C3*np.exp(x1) + C4*np.exp(x2);
        Tu2[i] = C1*np.exp(x3) + C2*np.exp(x4);
        Ty[i] = Td[i]
        Ty[nn-i - 1] = Tu[i]
        Ty2[i] = Td[i]
        Ty2[nn-i - 1] = Tu2[i]
        zy[i] = z[i]
        zy[nn-i-1] = z[i]
    return z,Td,Tu,zy,Ty,Q

def  Calcul_TBeir_coaxial(R12,R1,CCC,L,Rs,cas = 'annulus'):
    N12 = L/(CCC*R12)
    if cas == 'inner':
        Ns2 = L/(CCC*(R1 + Rs))
        Ns1 = 0
    elif cas == 'annulus':
        Ns1 = L/(CCC*(R1 + Rs))
        Ns2 = 0
    num = np.sqrt((Ns1 - Ns2)**2 + 4.0*((N12+Ns1)*(N12+Ns2) - N12**2))
    a1 = (-(Ns1-Ns2) + num)/2.0
    a2 = (-(Ns1-Ns2) - num)/2.0
    del1 = (Ns1 + N12 + a1)/N12
    del2 = (Ns1 + N12 + a2)/N12
    C1 = (del2-1)*np.exp(a2)/((del2-1)*np.exp(a2) - (del1- 1)*np.exp(a1))
    C2 = 1 - C1
    C3 = (Ns1 + N12 + a1)*C1/N12
    C4 = (Ns1 + N12 + a2)*C2/N12
    Q = 0
    dz = L/15.0
    z = np.arange(0,L+dz,dz)
    n = len(z)
    Td = np.zeros(n)
    Tu = np.zeros(n)
    Tu2 = np.zeros(n)
    nn = 2*n
    Ty = np.zeros(nn)
    Ty2 = np.zeros(nn)
    zy = np.zeros(nn)
    for i in range(0,n):
        x1 = a1*z[i]/L;
        x2 = a2*z[i]/L;
        x3 = a1*(2*L - z[i])/L
        x4 = a2*(2*L - z[i])/L
        Td[i] = C1*np.exp(x1) + C2*np.exp(x2)
        Tu[i] = C3*np.exp(x1) + C4*np.exp(x2)
        Tu2[i] = C1*np.exp(x3) + C2*np.exp(x4)
        Ty[i] = Td[i]
        Ty[nn-i - 1] = Tu[i]
        Ty2[i] = Td[i]
        Ty2[nn-i - 1] = Tu2[i]
        zy[i] = z[i]
        zy[nn-i-1] = z[i]
    return z,Td,Tu,zy,Ty,Q

def  Calcul_TBeir_coaxial_lin(R12,R1,CCC,L,Rs,ms,cas = 'annulus'):
    N12 = L/(CCC*R12)
    if cas == 'inner':
        Ns2 = L/(CCC*(R1 + Rs))
        Ns1 = 0
    elif cas == 'annulus':
        Ns1 = L/(CCC*(R1 + Rs))
        Ns2 = 0
    Delta = (N12+Ns1)*(N12+Ns2) - N12**2
    num = np.sqrt((Ns1 - Ns2)**2 + 4.0*((N12+Ns1)*(N12+Ns2) - N12**2))
    a1 = (-(Ns1-Ns2) + num)/2.0
    a2 = (-(Ns1-Ns2) - num)/2.0
    del1 = (Ns1 + N12 + a1)/N12
    del2 = (Ns1 + N12 + a2)/N12
    den = (del2-1)*np.exp(a2) - (del1- 1)*np.exp(a1)
    num1 = (del2-1)*np.exp(a2)
    num2 = (Ns1 + Ns2)*((a2*np.exp(a1) - a1*np.exp(a2))/(a1-a2) + 1 )/(a1*a2)
    C1s = num1/den
    C2s = 1 - C1s
    C3s = (Ns1 + N12 + a1)*C1s/N12
    C4s = (Ns1 + N12 + a2)*C2s/N12
    dens = 1 - (C3s + C4s)
    m = ms/(1 - ms*del1*num2/(den*dens))
    m = ms
    X = m*(Ns1 + Ns2)*((a2*np.exp(a1) - a1*np.exp(a2))/(a1-a2) + 1 )/(a1*a2)
    Y = m*(Ns1 + Ns2)/(a1*a2)*((a2*np.exp(a1) - a1*np.exp(a2))/(a1-a2) + 1)
    C1 = (num1-m*num2)/den
    C2 = 1 - C1
    C3 = (Ns1 + N12 + a1)*C1/N12
    C4 = (Ns1 + N12 + a2)*C2/N12
    Q = 0
    dz = L/15.0
    z = np.arange(0,L+dz,dz)
    n = len(z)
    Td = np.zeros(n)
    Tu = np.zeros(n)
    nn = 2*n
    Ty = np.zeros(nn)
    zy = np.zeros(nn)
    for i in range(0,n):
        x1 = a1*z[i]/L
        x2 = a2*z[i]/L
        W1 = (1 - np.exp(-x1))/a1
        W2 = (1 - np.exp(-x2))/a2
        Y1 = Delta*(1 - (1+x1)*np.exp(-x1))/a1**2
        Y2 = Delta*(1 - (1+x2)*np.exp(-x2))/a2**2
        V1 = Ns1*W1 - Y1
        V2 = Ns1*W2 - Y2
        V3 = -Ns2*W1 - Y1
        V4 = -Ns2*W2 - Y2
        I11 = np.exp(x1)*m*V1/(a1 - a2)
        I12 = np.exp(x2)*m*V2/(a2 - a1)
        I21 = np.exp(x1)*m*V3/(a1 - a2)
        I22 = np.exp(x2)*m*V4/(a2 - a1)
        Ix = Ns1*(np.exp(x1)/a1 - np.exp(x2)/a2 + 1/a2 - 1/a1)
        Iz = -Ns2*(np.exp(x1)/a1 - np.exp(x2)/a2 + 1/a2 - 1/a1)
        Iy = Delta*(np.exp(x1)/a1**2 - np.exp(x2)/a2**2 + (1+ x2)/a2**2 - (1+ x1)/a1**2)
        J = m/(a1-a2)*(Ix - Iy)
        K = m/(a1-a2)*(Iz - Iy)
#        Td[i] = C1*np.exp(x1) + C2*np.exp(x2) + I11 + I12
#        Tu[i] = C3*np.exp(x1) + C4*np.exp(x2) + I21 + I22
        Td[i] = C1*np.exp(x1) + C2*np.exp(x2) + J
        Tu[i] = C3*np.exp(x1) + C4*np.exp(x2) + K
        Ty[i] = Td[i]
        Ty[nn-i - 1] = Tu[i]
        zy[i] = z[i]
        zy[nn-i-1] = z[i]
    return z,Td,Tu,zy,Ty,Q




def Calcul_RaRb(Rb1,Ra1,x1,y1,CCC,L):


    def fct(x):
        y = zeros(2)
        Rb = x[0]
        Ra = x[1]
        eta = L/(CCC*np.sqrt(Ra*Rb))
        R12 = 4*Ra*Rb/(4*Rb-Ra)
        xsi = L/(CCC*2*Rb*eta)
        zeta = xsi*2*Rb/R12
        y[0] = (np.cosh(eta)-(2*xsi*np.sinh(eta)*zeta/(np.cosh(eta)+xsi*np.sinh(eta))+xsi)*np.sinh(eta)-x1)
        y[1] = ((np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta)+xsi*np.sinh(eta))-y1)
        return y
    xi = [Rb1,Ra1]
    y = fsolve(fct,xi)
    z = fct(y)
    Rb = y[0]
    Ra = y[1]
    return Rb,Ra,z

def Calcul_RaRb_Beirm(R1i,R12i,Rs,x1,y1,CCC,L):


    def fct(x):
        y = zeros(2)
        R1 = x[0]
        R12 = x[1]
        N12 = L/(CCC*R12)
        Ns1 = L/(CCC*(R1 + Rs))
        Ns2 = Ns1
        num = np.sqrt((Ns1 - Ns2)**2 + 4.0*((N12+Ns1)*(N12+Ns2) - N12**2))
        a1 = (-(Ns1-Ns2) + num)/2.0
        a2 = (-(Ns1-Ns2) - num)/2.0
        C1 = (Ns1 + a2)*exp(a2)/((Ns1+a2)*exp(a2) - (Ns1 + a1)*exp(a1))
        C2 = 1 - C1
        C3 = (Ns1 + N12 + a1)*C1/N12
        C4 = (Ns1 + N12 + a2)*C2/N12
        zz = C1*np.exp(a1)  + C2*np.exp(a2) - x1
        ww = C3 + C4-y1
        y = np.mean(zz)**2 + np.mean(ww)**2
        return y
    xi = [R1i,R12i]
    yy = minimize(fct,xi,method='Nelder-Mead')
    y = yy.x
    R1 = y[0]
    R12 = y[1]
    return R1,R12

def Calcul_RaRb_Beir(R1i,R12i,Rs,x1,y1,CCC,L):


    def fct(x):
        y = zeros(2)
        R1 = x[0]
        R12 = x[1]
        N12 = L/(CCC*R12)
        Ns1 = L/(CCC*(R1 + Rs))
        Ns2 = Ns1
        num = np.sqrt((Ns1 - Ns2)**2 + 4.0*((N12+Ns1)*(N12+Ns2) - N12**2))
        a1 = (-(Ns1-Ns2) + num)/2.0
        a2 = (-(Ns1-Ns2) - num)/2.0
        C1 = (Ns1 + a2)*exp(a2)/((Ns1+a2)*exp(a2) - (Ns1 + a1)*exp(a1))
        C2 = 1 - C1
        C3 = (Ns1 + N12 + a1)*C1/N12
        C4 = (Ns1 + N12 + a2)*C2/N12
        y[0] = C1*np.exp(a1)  + C2*np.exp(a2) - x1
        y[1] = C3 + C4-y1
        return y

    xi = [R1i,R12i]
    y = fsolve(fct,xi)
    R1 = y[0]
    R12 = y[1]
    return R1,R12


def Calcul_RaRbq(Rb1,Ra1,x1,x2,CCC1,CCC2,L):


    def fct(x):
        y = zeros(2)
        Rb = x[0]
        Ra = x[1]
        eta1 = L/(CCC1*sqrt(Ra*Rb))
        Rbs1 = Rb*eta1/tanh(eta1)
        eta2 = L/(CCC2*sqrt(Ra*Rb))
        Rbs2 = Rb*eta2/tanh(eta2)
        y = (Rbs1 - x1)**2 + (Rbs2 - x2)**2
        return y
    xi = [Rb1,Ra1]
    yy = minimize(fct,xi,method='Nelder-Mead')
    y = yy.x
    z = fct(y)
    Rb = y[0]
    Ra = y[1]
    return Rb,Ra,z

def Calcul_RaRbqq(Rb1,Ra1,x1,x2,CCC1,CCC2,L):


    def fct(x):
        y = zeros(2)
        Rb = x[0]
        Ra = x[1]
        eta1 = L/(CCC1*sqrt(Ra*Rb))
        Rbs1 = Rb*(1 + eta1*eta1/3.0)
        eta2 = L/(CCC2*sqrt(Ra*Rb))
        Rbs2 = Rb*(1 + eta2*eta2/3.0)
        y = (Rbs1 - x1)**2 + (Rbs2 - x2)**2
        return y
    xi = [Rb1,Ra1]
    yy = minimize(fct,xi,method='Nelder-Mead')
    y = yy.x
    z = fct(y)
    Rb = y[0]
    Ra = y[1]
    return Rb,Ra,z



def G_function_mils(Fo,Pe,rt = 1,the= 0):
    u = rt**2/(4*Fo)
    def g_integrand(x):
        b = rt*Pe
        y = np.exp(-x -b*b/(16*x))/x
        return y
    I = quad(g_integrand, u,Inf)
    G = I[0]*np.exp(rt*Pe*np.cos(the)/2)/(4*pi)
    return G

def G_function_mils_mean(Fo,Pe):
    u = 1/(4*Fo)
    def g_integrand(x):
        b = Pe
        y = np.exp(-x -b*b/(16*x))/x
        return y
    I = quad(g_integrand, u,Inf)
    G = I[0]*i0(Pe/2)/(4*pi)
    return G

def G_function_mfls(Fo,Pe,zt = 500,ht = 1000,rt=1,the =0,zob = 0):
    a = zob
    b = (Pe*rt/4)**2
    def funx(u):
        x = np.sqrt(u)/rt
        y = np.exp(-u-b/u)*(erf(x*(zt+a))+ erf(x*(zt-a)) -erf(x*(zt+ht+a))-erf(x*(zt-ht-a)))/u
        return y
    xi = rt**2/(4*Fo)
    y2 = quad(funx,xi,Inf)
    G = y2[0]*np.exp(rt*Pe*np.cos(the)/2)/(8*pi)
    return G

def G_function_mfls_mean(Fo,Pe,zt = 500,ht = 1000,zob = 0):
    a = zob
    b = (Pe/4)**2
    def funx(u):
        x = np.sqrt(u)
        y = np.exp(-u-b/u)*(erf(x*(zt+a))+ erf(x*(zt-a)) -erf(x*(zt+ht+a))-erf(x*(zt-ht-a)))/u
        return y
    xi = 1/(4*Fo)
    y2 = quad(funx,xi,Inf)
    G = y2[0]**i0(Pe/2)/(8*pi)
    return G

def G_function_mfls_bar(Fo,Pe,ht = 1000,zob = 0):
    a = zob
    b = (Pe/4)**2
    def funx(u):
        x = np.sqrt(u)
        g = 2*ierf(x*(ht+2*a)) + 2*ierf(x*ht) - ierf(x*2*a) - ierf(x*2*(ht+a))
        y = np.exp(-u-b/u)*g/u**(1.5)
        return y
    xi = 1/(4*Fo)
    y2 = quad(funx,xi,Inf)
    G = y2[0]*i0(Pe/2)/(8*pi*ht)
    return G

def leaky_function(u,rho):

    tau = np.log(rho/(2*u))
    if tau > 100:
        tau = 100
    h_inf = k0(rho)
    expintrho = expn(1,rho)
    w = (expintrho-h_inf)/(expintrho-expn(1,rho/2))
    ex1 = expn(1,rho/2*np.exp(abs(tau)))
    ex2 = expn(1,rho*np.cosh(tau))
    I = h_inf - w*ex1+(w-1)*ex2
    h=h_inf+np.sign(tau)*I
    return h

def Tcoax_inner(zz,xsi,gam,eta):
    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz = np.exp(gam*zz/2)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz = np.exp(gam*zz/2)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    the = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    return Tdz,Tuz

def Tcoax_inner_lin(zz,xsi,gam,eta,g):
    k12n = 2*xsi/(eta*(1-xsi**2))
    aa = (1 + xsi**2)/(2*xsi)
    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz2a = np.exp(gam*zz)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tdz2n = g*(zz + k12n*(Tdz2a-1) - np.exp(gam*(zz-1))/(eta*deltan)*np.sinh(eta*zz))
    Tdz2 = Tdz2a + Tdz2n
    Tuz2a = np.exp(gam*zz)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz2n = g*(zz + k12n*(-np.exp(gam*(zz-1))/deltan*(aa*np.sinh(eta*zz) + np.cosh(eta*zz)) + Tuz2a))
    Tuz2 = Tuz2a + Tuz2n
    theo = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    theon = theo + g*k12n*(theo - np.exp(-gam)/deltan)
    return Tdz2,Tuz2,theon
def Tcoax_innerb(zz,xsi,gam,eta,gx):
    g = 2*gx/(2+gx)
    k12n = 2*xsi/(eta*(1-xsi**2))
    aa = (1 + xsi**2)/(2*xsi)
    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz2a = np.exp(gam*zz)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tdz2n = g*(zz + k12n*(Tdz2a-1) - np.exp(gam*(zz-1))/(eta*deltan)*np.sinh(eta*zz))
    Tdz2 = Tdz2a + Tdz2n
    x = gx/2
    Tdz2 = Tdz2*(1 + x) - x
    Tuz2a = np.exp(gam*zz)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz2n = g*(zz + k12n*(-np.exp(gam*(zz-1))/deltan*(aa*np.sinh(eta*zz) + np.cosh(eta*zz)) + Tuz2a))
    Tuz2 = Tuz2a + Tuz2n
    Tuz2 = Tuz2*(1 + x) - x
    theo = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    theon = theo + g*k12n*(theo - np.exp(-gam)/deltan)
    return Tdz2,Tuz2,theon

def Tcoax_annulus(zz,xsi,gam,eta):

    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz = np.exp(-gam*zz/2)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz = np.exp(-gam*zz/2)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    theo = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    return Tdz,Tuz
def Tcoax_annulus_lin(zz,xsi,gam,eta,g=0):

    k12n = 2*xsi/(eta*(1-xsi**2))
    aa = (1 + xsi**2)/(2*xsi)
    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz2a = np.exp(-gam*zz)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tdz2n = g*(zz - np.sinh(eta*zz)*np.exp(gam*(1-zz))/(eta*deltan))
    Tdz2 = Tdz2a + Tdz2n
    Tuz2a = np.exp(-gam*zz)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz2n = g*(zz + k12n - k12n/deltan*np.exp(gam*(1-zz))*(np.sinh(eta*zz)*aa + np.cosh(eta*zz)))
    Tuz2 = Tuz2a + Tuz2n
    theo = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    theon = theo + g*k12n - g*k12n*np.exp(gam)/deltan
    return Tdz2,Tuz2,theon

def Tcoax_annulusb(zz,xsi,gam,eta,gx):
    g = 2*gx/(2+gx)
    k12n = 2*xsi/(eta*(1-xsi**2))
    aa = (1 + xsi**2)/(2*xsi)
    deltan = (np.cosh(eta)+xsi*np.sinh(eta));
    Tdz2a = np.exp(-gam*zz)*(np.cosh(eta*(1-zz))+xsi*np.sinh(eta*(1-zz)))/deltan
    Tdz2n = g*(zz - np.sinh(eta*zz)*np.exp(gam*(1-zz))/(eta*deltan))
    Tdz2 = Tdz2a + Tdz2n
    x = gx/2.0
    Tdz2 = Tdz2*(1 + x) - x
    Tuz2a = np.exp(-gam*zz)*(np.cosh(eta*(1-zz))-xsi*np.sinh(eta*(1-zz)))/deltan
    Tuz2n = g*(zz + k12n - k12n/deltan*np.exp(gam*(1-zz))*(np.sinh(eta*zz)*aa + np.cosh(eta*zz)))
    Tuz2 = Tuz2a + Tuz2n
    Tuz2 = Tuz2*(1 + x) - x
    theo = (np.cosh(eta)-xsi*np.sinh(eta))/(np.cosh(eta) + xsi*np.sinh(eta))
    theon = theo + g*k12n - g*k12n*np.exp(gam)/deltan
    return Tdz2,Tuz2,theon

def Tcoax_inner_flux(zz,R1,R12):
    Tuz = R1 - 0.5 - 1/(6*R12) +  zz + 1/(2*R12)*(1-zz)**2
    Tdz = R1 + 0.5 - 1/(6*R12) + 1/(2*R12)*(1-zz)**2
    theo = R1 - 0.5  + 1/(3*R12)
    thei = R1 + 0.5  + 1/(3*R12)
    return Tdz,Tuz,theo,thei


def Tcoax_annulus_flux(zz,R1,R12):
    Tdz = R1 + 0.5 - 1/(6*R12) -  zz + 1/(2*R12)*(1-zz)**2
    Tuz = R1 - 0.5 - 1/(6*R12) + 1/(2*R12)*(1-zz)**2
    theo = R1 - 0.5  + 1/(3*R12)
    thei = R1 + 0.5  + 1/(3*R12)
    return Tdz,Tuz,theo,thei

def Tcoax_annulus_flux_lin(zz,R1,R12):
    Tdz = R1 - 1/3  - 1/(12*R12) + (1-  zz)**2  + 1/(3*R12)*(1-zz)**3
    Tuz = R1 - 1/3 - 1/(12*R12) + 1/(3*R12)*(1-zz)**3
    theo = R1 -  1/3  + 1/(4*R12)
    thei = R1 + 2/3  + 1/(4*R12)
    return Tdz,Tuz,theo,thei

def Tcoax_inner_flux_lin(zz,R1,R12):
    Tuz = R1 - 1/3 - 1/(4*R12) + zz**2  + (zz**3 - 3*zz+2)/(3*R12)
    Tdz = R1 + 2/3 - 1/(4*R12)  + (zz**3 - 3*zz+2)/(3*R12)
    theo = R1 -  1/3  + 5/(12*R12)
    thei = R1 + 2/3  + 5/(12*R12)
    return Tdz,Tuz,theo,thei


