#
#
from finance_md import *
import  numpy as np

# a
mib = 2200
ti = 0.07
inte = 0.04
gain = 800/(1 + ti) + 1000/(1+ti)**2  +  1200/(1+ti)**3
cout = 50/(1 + ti) + 50*(1 + inte)/(1+ti)**2  +  50*(1 + inte)**2/(1+ti)**3
vana = -mib + gain - cout
den = pwf(3,0,ti)
EAC1 = vana/den
print(EAC1)

gaina = future_value(800,2,ti) + future_value(1000,1,ti)  +  1200
cout = future_value(50,2,ti) + future_value(50*(1+inte),1,ti)  +  50*(1+inte)**2
mibf = future_value(mib,3,ti)
den = fwf(3,0,ti)
vfna = -mibf + gaina - cout
EAC2 = vfna/den
print(EAC2)
