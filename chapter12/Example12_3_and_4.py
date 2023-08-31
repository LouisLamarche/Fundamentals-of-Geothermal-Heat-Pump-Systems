#
#
from finance_md import *
import  numpy as np

mi = 10000.0        # cout total
Npret = 3
N = 3
t = 0.08    # taux d'actualisation
inte = 0.06
paiement = mi/pwf(Npret,0,inte)     # paiment annuel
inter  = np.zeros(N)
cap = np.zeros(N)
reste = mi
va1 = 0.0
for i in range(0,N):
    inter[i] = inte*reste
    cap[i] =  paiement - inter[i]
    reste = reste - cap[i]
    print(reste)
    va1 = va1 + present_value(inter[i],i+1,0,t)
va2 = mi*pwf_int(Npret,N,inte,t)
print ('The discounted value of the 2 interest are:  ',va1)
print ('The discounted value of the 2 interest are:  ',va2)
