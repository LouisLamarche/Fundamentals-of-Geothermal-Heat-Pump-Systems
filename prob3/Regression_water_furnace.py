from numpy import *
from conversion_md import *
from performance_water_furnace import *
from scipy.optimize import curve_fit,newton,minimize,fsolve
from  matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

Tlv = arange(10,36,2)
Tsv = arange(-1,20,2)
ni = len(Tlv)
nj = len(Tsv)
nt = ni*nj
HC = zeros(nt)
POW = zeros(nt)
CC = zeros(nt)
COP = zeros(nt)
x = zeros((2,nt))
X,Y = meshgrid(Tlv,Tsv)
Z2 = zeros((nj,ni))
Z3 = zeros((nj,ni))
Z4 = zeros((nj,ni))
k = 0
for i in range(0,ni):
    for j in range(0,nj):
        x[0,k] = Tlv[i]
        x[1,k] = Tsv[j]
        TL = F_C(Tlv[i])
        TS = F_C(Tsv[j])
        CC_cap,Pow = performance_water_furnace(TS,TL)
        CC[k] = W_BTUhr(CC_cap)
        Z2[j,i] = W_BTUhr(CC_cap)
        POW[k] = Pow
        COP[k] = W_BTUhr(CC_cap)/Pow
        HC[k] =  W_BTUhr(CC_cap) + Pow
        Z3[j,i] =  COP[k]
        Z4[j,i] = HC[k]
        k = k+1
def modele(x,a,b,c,d,e):
    y = a+ b*x[0] + c*x[1] + d*x[0]*x[1]  + e*x[1]*x[1]
    return y
p1a,resn = curve_fit(modele,x,CC)
p1 = [11.02, -0.033,-0.164, 0.012, -4.86e-5]
Y2 = p1[0]+ p1[1]*X+ p1[2]*Y + p1[3]*X*Y + p1[4]*Y*Y
p2a,resn = curve_fit(modele,x,HC)
p2 = [12.38, -0.04,-0.136, 0.012, -4.43e-5]
Y4 = p2[0]+ p2[1]*X+ p2[2]*Y   + p2[3]*X*Y + p2[4]*Y*Y
p3,resn = curve_fit(modele,x,COP)
Y3 = p3[0]+ p3[1]*X+ p3[2]*Y   + p3[3]*X*Y + p3[4]*Y*Y
fig = figure()
ax = fig.gca(projection='3d')
# Plot the surface.
surf = ax.plot_surface(X, Y, Y4, cmap=cm.coolwarm)
surf = ax.scatter(X, Y, Z4, cmap=cm.coolwarm)
show()
exit()
