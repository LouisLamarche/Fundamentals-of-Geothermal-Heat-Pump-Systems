#
#
from finance_md import *
import numpy as np
c1b = np.array([3240,3200,3150,3060,2900])
mi1 = 12000
c2b = np.array([3780,3790,3780,3740,3670])
mi2 = 10000
CV1= mi1 + sum(c1b)
CV2= mi2 + sum(c2b)
t = 0.08    # taux d'actualisation
ccv1 = mi1
ccv2 = mi2
n = len(c1b)
pv1 = np.zeros(n)
pv2 = np.zeros(n)
for i in range(0,n):
    pv1[i] = c1b[i]/(1+t)**(i+1)
    ccv1 = ccv1 + c1b[i]/(1+t)**(i+1)
    pv2[i] = c2b[i]/(1+t)**(i+1)
    ccv2 = ccv2 + c2b[i]/(1+t)**(i+1)
print('LCC project 1 = ',ccv1)
print('LCC project 2 = ',ccv2)
