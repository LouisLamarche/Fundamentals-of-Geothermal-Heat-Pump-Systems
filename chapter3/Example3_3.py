#coding: utf-8
# exemple 3.3
# written by Louis Lamarche
# 15 august 2017
#
import  numpy as np
from conversion_md import *
from CoolProp.CoolProp import *
from tabulate import tabulate

# data
def KC(x):
    return x - 273.15
def CK(x):
    return x + 273.15
#

# data

fluid = 'R134a'
patm = 101.325*1000.0
p1 = 1.70e5 # pascal
p2 = 1000.0*1000 # pascal
p3 = 995.0*1000 # pascal
p4 = 1.80e5 # pascal
p5 = 1.75e5 # pascal
etas = 0.7
#
T4k = PropsSI('T', 'P', p4, 'Q', 1, fluid)
T4c = KC(T4k)
T5ks = PropsSI('T', 'P', p5, 'Q', 1, fluid)
T5s = KC(T5ks)
T5k = T5ks + 4.0
T5c = KC(T5k)
h5 = PropsSI('H', 'T',T5k,'P', p5, fluid)
s5 = PropsSI('S', 'T',T5k,'P', p5, fluid)
T1k = T5k + 1.0
T1c = T5c + 1.0
h1 = PropsSI('H', 'T',T1k,'P', p1, fluid)
s1 = PropsSI('S', 'T',T1k,'P', p1, fluid)
h2s = PropsSI('H', 'P', p2, 'S', s1, fluid)
h2 = h1 +  (h2s - h1)/etas
s2 = PropsSI('S', 'P', p2, 'H', h2, fluid)
T2k = PropsSI('T', 'P', p2, 'H', h2, fluid)
T2c = KC(T2k)
T3ks = PropsSI('T', 'P', p3, 'Q', 0, fluid)
T3s = KC(T3ks)
T3k = T3ks - 3
T3c = KC(T3k)
h3 = PropsSI('H', 'T',T3k,'P', p3, fluid)
s3 = PropsSI('S', 'T',T3k,'P', p3, fluid)
h4 = h3
x4 = PropsSI('Q', 'P', p4, 'H', h4, fluid)
s4 = PropsSI('S', 'P', p4, 'Q', x4, fluid)
qin = h5 - h4
qout = h2 - h3
win = h2  - h1
COP =qin/win
Qin = 5000.0
mp = Qin/qin
print ('COP = ', '%.2f' % COP)
print ('mdot = ','%.3f' % mp, ' kg/s')
head = ['point','T(C)', 'p(kPa)','h(kJ/kg)','s(kJ/kg-K)']
pt1 = ['1','%.1f' % T1c, '%.1f' % (p1/1000),  '%.0f' % (h1/1000), '%.3f' % (s1/1000)]
pt2 = ['2','%.1f' % T2c, '%.1f' % (p2/1000),  '%.0f' % (h2/1000), '%.3f' % (s2/1000)]
pt3 = ['3','%.1f' % T3c, '%.1f' % (p3/1000),  '%.0f' % (h3/1000), '%.3f' % (s3/1000)]
pt4 = ['4','%.1f' % T4c, '%.1f' % (p4/1000),  '%.0f' % (h4/1000), '%.3f' % (s4/1000)]
pt5 = ['5','%.1f' % T5c, '%.1f' % (p5/1000),  '%.0f' % (h5/1000), '%.3f' % (s5/1000)]
print(tabulate([pt1,pt2,pt3,pt4,pt5], headers= head))
