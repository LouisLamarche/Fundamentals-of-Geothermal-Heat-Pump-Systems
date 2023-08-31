#coding: utf-8
#
import numpy as np
from conversion_mod import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
from tabulate import tabulate

# data
def KC(x):
    return x - 273.15
def CK(x):
    return x + 273.15
#
# data
#
fluid = 'R134a'
patm = 101.325*1000.0
p1 = p4 = 1.8e5
p2 = p3 = 1e6
etac = 0.8
surch = 5
#
#
T4k = PropsSI('T', 'P', p4, 'Q', 1, fluid)
T4c = KC(T4k)
T1c = T4c + surch
T1k = T4k + surch
h1 = PropsSI('H', 'P', p1, 'T', T1k, fluid)
h1 = PropsSI('H', 'P', p1, 'T', T1k, fluid)
s1 = PropsSI('S', 'P', p1, 'T', T1k, fluid)
s2s = s1
h2s = PropsSI('H', 'P', p2, 'S', s2s, fluid)
h2 = h1 + (h2s - h1)/etac
T2k = PropsSI('T', 'P', p2, 'H', h2, fluid)
s2 = PropsSI('S', 'P', p2, 'H', h2, fluid)
T2c = KC(T2k)
T3k = PropsSI('T', 'P', p3, 'Q', 0, fluid)
T3c = KC(T3k)
h3 = PropsSI('H', 'P', p3, 'Q', 0, fluid)
s3 = PropsSI('S', 'P', p3, 'Q', 0, fluid)
h4 = h3
x4 = PropsSI('Q', 'P', p4, 'H', h4, fluid)
s4 = PropsSI('S', 'P', p1, 'Q', x4, fluid)
qin = h1 - h4
qout = h2 - h3
win = h2  - h1
COPC =qin/win
COPH =qout/win
print ('COP (cooling) = ', '%.2f' % COPC)
print ('COP (heating) = ', '%.2f' % COPH)
head = ['point','T(C)', 'p(kPa)','h(kJ/kg)','s(kJ/kg-K)']
pt1 = ['1','%.1f' % T1c, '%.1f' % (p1/1000),  '%.0f' % (h1/1000), '%.2f' % (s1/1000)]
pt2 = ['2','%.1f' % T2c, '%.1f' % (p2/1000),  '%.0f' % (h2/1000), '%.2f' % (s2/1000)]
pt3 = ['3','%.1f' % T3c, '%.1f' % (p3/1000),  '%.0f' % (h3/1000), '%.2f' % (s3/1000)]
pt4 = ['4','%.1f' % T4c, '%.1f' % (p4/1000),  '%.0f' % (h4/1000), '%.2f' % (s4/1000)]
print(tabulate([pt1,pt2,pt3,pt4], headers= head))
