#
#
from finance_md import *
import numpy as np
from scipy.optimize import newton

# a
inte = 0.04
gainv = np.array([800,1000,1200])
costs = np.array([50,50,50])
costv = np.array([50,50*(1+inte),50*(1+inte)**2])
profits = gainv - costv
iir = 0.13007979152900095
mi = 2200
G = sum(gainv)
C = 50*fwf(3,0,inte)
P1 = mi*(1+iir)**3 - mi
print('profit 1= ',P1)
P2 = G - C - mi
print('profit 2= ',P2)
P3a = - mi - 50*fwf(3,inte,iir) + gainv[0]*(1 + iir)**2 + gainv[1]*(1 + iir) +gainv[2]
P3b = - mi +  profits[0]*(1 + iir)**2 + profits[1]*(1 + iir) + profits[2]
g = (mi + costv[0]/(1+iir) + costv[1]/(1+iir)**2 + costv[2]/(1+iir)**3)*(1+iir)**3
d = gainv[0]*(1 + iir)**2 + gainv[1]*(1 + iir) +gainv[2]
print('profit 3 = ',P3a,P3b)
