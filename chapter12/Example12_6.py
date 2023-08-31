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
n = len(cbb)
gaind1 = np.zeros(n)
npv = -mi1
gains1 = cbb - c1b
for i in range(0,n):
    gaind1[i] = gains1[i]/(1+t)**(i+1)
    npv = npv + gaind1[i]
print('Net present value = ',npv)

