#
#
from finance_md import *
import numpy as np
from matplotlib.pyplot import *
cbb = np.array([5400,5540,5670,5780,5880])
c1b = np.array([3240,3200,3150,3060,2900])
mi1 = 12000
CVb = np.sum(cbb)
CV1= mi1 + np.sum(c1b)
t = 0.08    # taux d'actualisation
ccvb = 0
ccv1 = mi1
n = len(cbb)
pvb = np.zeros(n)
pv1 = np.zeros(n)
gaind1 = np.zeros(n)
cub = np.zeros(n+1)
cu1 = np.zeros(n+1)
cnb = np.zeros(n+1)
cn1 = np.zeros(n+1)
cu1[0] = ccv1
cn1[0] = ccv1
for i in range(0,n):
    pvb[i] = cbb[i]/(1+t)**(i+1)
    ccvb = ccvb + cbb[i]/(1+t)**(i+1)
    cub[i+1] = ccvb
    cnb[i+1] = cnb[i] + cbb[i]
    pv1[i] = c1b[i]/(1+t)**(i+1)
    ccv1 = ccv1 + c1b[i]/(1+t)**(i+1)
    cu1[i+1] = ccv1
    cn1[i+1] = cn1[i] + c1b[i]
print(CVb,CV1)
print(ccvb,ccv1)
gains1 = ccv1 - ccvb
print('Cost of new project = ',gains1)

