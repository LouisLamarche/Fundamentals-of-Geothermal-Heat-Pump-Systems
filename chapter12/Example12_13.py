#
#
from finance_md import *
import numpy as np
from scipy.optimize import newton

# a
inte = 0.04
gainv = np.array([800,1000,1200])
costv = np.array([50,50*(1+inte),50*(1+inte)**2])
flow = gainv - costv
profits = np.maximum(gainv - costv,0)
pertes = np.maximum(costv - gainv,0)
iir = 0.13007979152900095
dpa = 2200
fr = 0.05      # finacial rate
rr = 0.07   # reinvestment rate
#fr = 0.13007979152900095
#rr = 0.13007979152900095
num= 0
n = len(profits)
m = len(pertes)
for i in range(0,n):
    num = num + profits[i]*(1 + rr)**(n-1-i)
den = dpa
for i in range(0,m):
    den = den + pertes[i]/(1 + fr)**(i+1)
miira1 = (num/den)**(1/3) - 1
miira2 = mirr(dpa,profits,rr,fr)
print ('miir is = ',miira1,miira2)
# b)
miirb = mirr(dpa,flow,0,fr)
print ('miir is = ',miirb)