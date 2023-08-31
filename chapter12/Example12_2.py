#
#
from finance_md import *
import  numpy as np

mi = 10000.0        # cout total
r = 0.10
mf1 = mi*(1+r)
print ('INterest compounded at the end of the year  ',mf1)
mfi = mi*(1+r/2)
mf2 = mfi*(1+r/2)
print ('INterest compounded twice a year  ',mf2)
N = 2
C = 2
req = (1 + r/N)**C - 1
mf3 = mi*(1+req)
print ('Equivalent interest rate  ',req)
print ('Interest compounded twice a year  ',mf3)
