#coding: latin-1
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse


def fct(y,n):
    x = (y-6+n)/(2*n -6)
    return x

y = 0.5
x1 = fct(y*6,4)
x2 = fct(y*6,5)

X = 0.4
Z = fct_lemire(1,X)
p1 = fct_lemire(1,X)/Z*y
p2 = fct_lemire(y,X)/Z
p3 = fct_lemire(4/6,X)/Z*x1 + fct_lemire(2/6,X)/Z*(1-x1)
p4 = fct_lemire(5/6,X)/Z*x2 + fct_lemire(1/6,X)/Z*(1-x2)
print(p1,p2,p3,p4)
print((p1+p2+p3+p4)/4)
#
# b
#
p1b = 0.5
p2b = y*((1-X)*y**2 + X)
y1 = 4/6
y2 = 2/6
p3b = x1*y1*((1-X)*y1**2 + X) + x2*y2*((1-X)*y2**2 + X)
y1 = 5/6
y2 = 1/6
p4b = x1*y1*((1-X)*y1**2 + X) + x2*y2*((1-X)*y2**2 + X)
etam1 = motor_efficiency(1.5,1,case = 'H')
eta_v1 = VFD_efficiency(3,1)
eta = etam1*eta_v1
etam2 = motor_efficiency(1.5,y,case = 'H')
eta_v2 = VFD_efficiency(3,y)
etapl = etam2*eta_v2
#p2b = p2b*eta/etapl
print(p1b,p2b,p3b,p4b)
print((p1b+p2b+p3b+p4b)/4)




