#
#
from finance_md import *
import  numpy as np

# a
mib = 2200
ti = 0.07
inte = 0.04
gaina = future_value(800,2,ti) + future_value(1000,1,ti)  +  1200
cout = future_value(50,2,ti) + future_value(50*(1+inte),1,ti)  +  50*(1+inte)**2
mibf = future_value(mib,3,ti)
vfna = -mibf + gaina - cout
print(vfna)
gainb = future_value(1000,2,ti) + future_value(1000,1,ti)  +  1000
vfnb = -mibf + gainb - cout
vfnb = -mibf + 1000*fwf(3,0,ti) - 50*fwf(3,inte,ti)
print(vfnb)
