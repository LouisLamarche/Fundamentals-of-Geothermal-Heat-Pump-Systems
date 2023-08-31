#
#
from finance_md import *
import numpy as np
from scipy.optimize import newton

# a
mi = 2200
inte = 0.04
def calcul_iir(t):
    gain = 800/(1 + t) + 1000/(1+t)**2  +  1200/(1+t)**3
    cout = 50/(1 + t) + 50*(1 + inte)/(1+t)**2  +  50*(1 + inte)**2/(1+t)**3
    npv = -mi + gain - cout
    return npv
ti = newton(calcul_iir,0.7)
# verification
print ('Thi iir (a) is  = ',ti)
gain = 800/(1 + ti) + 1000/(1+ti)**2  +  1200/(1+ti)**3
cout = 50/(1 + ti) + 50*(1 + inte)/(1+ti)**2  +  50*(1 + inte)**2/(1+ti)**3
npv = -mi + gain - cout
print ('The npv with t = iir =  ',npv)

