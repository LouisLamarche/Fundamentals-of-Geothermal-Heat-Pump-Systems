#
#
from finance_md import *
import numpy as np

# a
mib = 2200
ti = 0.07
inte = 0.04
gain = 800/(1 + ti) + 1000/(1+ti)**2  +  1200/(1+ti)**3
cout = 50/(1 + ti) + 50*(1 + inte)/(1+ti)**2  +  50*(1 + inte)**2/(1+ti)**3
vana = -mib + gain - cout
print(vana)
gain = 1000/(1 + ti) + 1000/(1+ti)**2  +  1000/(1+ti)**3
cout = 50/(1 + ti) + 50*(1 + inte)/(1+ti)**2  +  50*(1 + inte)**2/(1+ti)**3
vanb = -mib + gain - cout
vanb = -mib + 1000*pwf(3,0,ti) - 50*pwf(3,inte,ti)
print(vanb)
