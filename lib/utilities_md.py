import numpy as np
import  scipy.special as sp
pi = np.pi
delta = np.finfo(np.float64).eps
def T_function(p,z,beta):
    z = np.sqrt(p)
    rho = sp.k0(z)/(z*spk1(z))
    y = np.exp(-z/(rho + beta))/p
    return y

def lgwt(N,a,b):
    N=N-1
    N1=N+1
    N2=N+2
    xu= np.linspace(-1,1,N1)
    ## Initial guess
    xi = np.arange(0,N1)
    y = np.cos((2*xi+1)*pi/(2*N+2))+(0.27/N1)*np.sin(pi*xu*N/N2)
    # Legendre-Gauss Vandermonde Matrix
    L = np.zeros((N1,N2))
    # Derivative of LGVM
    Lp = np.zeros((N1,N2))
    fl = np.ones(N1)
    L[:,0] = fl
    Lp[:,0] = fl

    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method
    y0 = 2
    # Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0)) > delta:
        L[:,1] =y
        for k in range(2,N2):
            L[:,k] = ((2*k-1)*y*L[:,k-1]-(k-1)*L[:,k-2])/k
        Lp = (N2)*(L[:,N1-1]-y*L[:,N2-1])/(1-y**2)
        y0=y
        y=y0 - L[:,N2-1]/Lp
    # Linear map from[-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2
    w = (b-a)/((1-y**2)*Lp**2)*(N2/N1)**2
    return x,w

def  cherche_index(xi,x):
    """ cherche l'index où x(i) <= xi < x(i+1)"""
    err = 0
    ok = False
    i = 0
    if (xi < x[0] or xi > x[len(x)-1]):
        err =-1
        i = np.nan
    elif xi == x[len(x)-1]:
        i = len(x)-2
    else:
        while not ok:
            if (xi>=x[i]) & (xi<x[i+1]):
                ok = True
            else:
                i=i+1
                if i >= len(x):
                    ok = True
                    err =-1
                    i = nan
    return i


def  cherche_index2(xi,x):
    """ cherche l'index où x(i) <= xi < x(i+1)"""
    err = 0
    ok = 1
    i = 0
    if (xi>0 and xi <= x[0] ):
        i = 0
    elif xi == x[len(x)-1]:
        i = len(x)-1
    else:
        while ok:
            if (xi>=x[i]) & (xi<x[i+1]):
                ok = 0
            else:
                i=i+1
                if i >= len(x):
                    ok =0
                    err =-1
                    i = nan
    return i
def mon_interp2(x,y,z,x1,y1):

    i = cherche_index(x1,x)
    j = cherche_index(y1,y)
    if isinstance(i,int) and isinstance(i,int):
        A = np.array([[x[i],y[j],x[i]*y[j],1],\
        [x[i+1],y[j],x[i+1]*y[j],1],\
        [x[i],y[j+1],x[i]*y[j+1],1],\
        [x[i+1],y[j+1],x[i+1]*y[j+1],1]])
        B = np.array([z[i,j],z[i+1,j],z[i,j+1],z[i+1,j+1]]);
        d = np.linalg.solve(A,B)
        z1 = d[0]*x1+d[1]*y1+d[2]*x1*y1+d[3]
        return z1
    else:
        print('No extrapolation')
        return -99

def combinaisons(n,k):

    if k == 0:
        C = 1
    elif k ==1:
        C = n
    else:
        p = n-k+1
        for i in range (2,k+1):
            p = p*(n-k+i)/(i*1.0)
        C = p
    return C

def  fac_log(n):
    if n==0:
        p=0
    else:
        p = np.log(n)
        i = n-1
        while (i!=0):
            p = p+np.log(i)
            i = i-1
    return np.longdouble(p)

def der_sec(x,f,ai=0,af=0,it=[0,0]):
    n = len(x)
    A = np.zeros((n,n))
    d = np.zeros(n)
    if it[0]==1:
       df=f[1]-f[0]
       dx=x[1]-x[0]
       d[0]=3*(df/dx-ai)/dx
       A[0,0]=1.0
       A[0,1]=0.5
    else:
       d[0]=ai
       A[0,0]=1
    if it[1]==1:
       df=f[n-1]-f[n-2]
       dx=x[n-1]-x[n-2]
       d[n-1]=3.0*(af-df/dx)/dx
       A[n-1,n-1]=1.0
       A[n-1,n-2] = 0.5
    else:
       d[n-1]=ai
       A[n-1,n-1]=1.0

    for i in range(1,n-1):
       A[i,i-1]=(x[i]-x[i-1])/6.0
       A[i,i]=(x[i+1]-x[i-1])/3.0
       A[i,i+1]=(x[i+1]-x[i])/6.0
       d[i]=(f[i+1]-f[i])/(x[i+1]-x[i])-(f[i]-f[i-1])/(x[i]-x[i-1])
    g = np.linalg.solve(A,d)
    return g

def  f_spline(x2,g,x,f):
    from mmath_mod import cherche_index
    n1 = len(x)-1
    xi = x2
    if  (xi <= x[n1]) and (xi>=x[0]):   # interpolation
        i = cherche_index(xi,x)
        d1  = xi-x[i]
        d2  = x[i+1]-xi
        dx = x[i+1]-x[i]
        y = (g[i]*d2**3+g[i+1]*d1**3)/(6.0*dx)+(f[i]/dx-dx*g[i]/6.0)*d2+(f[i+1]/dx-dx*g[i+1]/6.0)*d1
    elif xi < x[0]:                       # extrapolation
        d1 = f[1] - f[0]
        d2 = x[1] - x[0]
        dx = xi-x[0]
        y = f[0]+(d1/d2-g[0]*d2/3.0+g[1]*d2/3.0)*dx+g[0]*dx*dx/2.0
    else:
        d1 = f[n1]-f[n1-1]
        d2 = x[n1]-x[n1-1]
        dx = xi-x[n1]
        y = f[n1]+(d1/d2-g[n1]*d2/3.0+g[n1-1]*d2/3.0)*dx+g[n1]*dx*dx/2.0
    return y

class cls_inv_laplace:


    def __init__(self,L):
        self.v =  np.longdouble(np.zeros(L))


    def gavsteh_log_spline(self,t,L,flag,ders,xx,gg):

        if flag == 0:
            nn2 = L//2
            for n in range(0,L):
                nn = n+1
                z = 0.0
                for k in range(int(np.floor((nn+1)/2.0)),min(nn,nn2)+1):
                    dz = nn2*np.log(k)+fac_log(2*k)-fac_log(k)-fac_log(nn2-k)-fac_log(k-1)-fac_log(nn-k)-fac_log(2*k-nn)
                    z = z + np.exp(dz)
                self.v[n] = (-1)**(nn+nn2)*z
        sum = 0.0
        ln2_on_t = np.log(2.0)/t
        for n in range (0,L):
            p = (n+1)*ln2_on_t
            pl = np.log(p)
            sum = sum + self.v[n]*f_spline(pl,ders,xx,gg)
        ilt = sum*ln2_on_t
        return ilt


    def gavsteh(self,t,L,flag,fct,**kwargs):

        if flag == 0:
            nn2 = L//2
            for n in range(0,L):
                nn = n+1
                z = 0.0
                for k in range(int(np.floor((nn+1)/2.0)),min(nn,nn2)+1):
                    dz = nn2*np.log(k)+fac_log(2*k)-fac_log(k)-fac_log(nn2-k)-fac_log(k-1)-fac_log(nn-k)-fac_log(2*k-nn)
                    z = z + np.exp(dz)
                self.v[n] = (-1)**(nn+nn2)*z
        sum = 0.0
        ln2_on_t = np.log(2.0)/t
        for n in range (0,L):
            p = (n+1)*ln2_on_t
            sum = sum + self.v[n]*fct(p,**kwargs)
        ilt = sum*ln2_on_t
        return ilt

