#
#
# Example 1.6 Kalina cycle
"""
 To use this program, a valid refprop dll should be installed on the PC
 By default it is at :
    c:\\Program Files (x86)\\REFPROP
 if not please change the module refprop_md
 """
#
import numpy as np
from CoolProp.CoolProp import  *
from scipy.optimize import newton
from refprop_md import *

T2p = 160.0

mp = 12     # refrigerant flow rate
mpw  = 80.0   # water flow rate
effe = 0.80
y = 0.8
pH = 3e6    # high pressure
p5 = p6 = p7 = p1 = p8 = p9 = pH
y3 = y4 = y5 = y6 = y7  = y
#
#
# debut

#
# hypothese point
#
T7 = 150
T7k = T7 + 273.15
a = calcul_state(T7k,p7,y7)
x7 = calcul_q_Tpx(T7k,p7,y7)
f = 1 - x7   #  mass fraction in the recuperator

#
# point 4 saturated liquid
#
T4 = 40
T4k = T4 + 273.15
p4 = PropsSI_nh3_h20('P','T',T4k,'Q',0,y4)
pL = p4      # low pressure
p2 = p3 = p10 = pL
h4 = PropsSI_nh3_h20('H','T',T4k,'Q',0,y4)
s4 = PropsSI_nh3_h20('S','T',T4k,'Q',0,y4)
#
# point 5 isentropic pump
#
s5 = s4
T5k = PropsSI_nh3_h20('T','P',p5,'S',s5,y5)
T5 = T5k - 273.15
h5 = PropsSI_nh3_h20('H','T',T5k,'S',s5,y5)
#
# point 7 inlet of separator
#
h7 = PropsSI_nh3_h20('H','P',p7,'Q',x7,y7)
T8k = T1k = T7k
T8 = T1 = T7
#
# point 1 dew point
y1 = calcul_x_Tpd(T7k,p7)
#
#y1 = 1.0  # Refprop a de probleme si y pres de 1.0
# on met y1 à 1 car refprop fonctionne mal pour x pres de 0.99
h1 = PropsSI_nh3_h20('H','T',T1k,'Q',1,y1)
s1 = PropsSI_nh3_h20('S','T',T1k,'Q',1,y1)
#
# point 8 bubble point
#
y8 = calcul_x_Tpb(T8k,p8)
fb = (y7 - y1)/(y8 - y1)
y9 = y10 = y8
h8 = (h7 - x7*h1)/f


# point 2s
s2 = s1
y2 = y1
h2s = PropsSI_nh3_h20('H','P',p2,'S',s2,y2)
etat = 1.0
h2 = h1 - (h1-h2s)*etat
T2k = PropsSI_nh3_h20('T','P',p2,'H',h2,y2)
T2 = T2k - 273.15

#
# point 9 et 6 sortie du récupérateur
#
h6m = PropsSI_nh3_h20('H', 'T',T8k,'P',p6,y6)
h9m = PropsSI_nh3_h20('H', 'T',T5k,'P',p9,y9)
qmax = min((h6m - h5),f*(h8 - h9m))
qech = effe*qmax
h6 = h5 + qech
h9 = h8 - qech/f
#
T6k = PropsSI_nh3_h20('T','P',p6,'H',h6,y6)
etat6 = calcul_state(T6k,p6,y6)
if etat6 == 'mix':
    x6 = calcul_q_Tpx(T6k,p6,y6)
else:
    x6 = 0
T6 = T6k - 273.15
T9k = PropsSI_nh3_h20('T','P',p9,'H',h9,y9)
T9 = T9k - 273.15
#
#
# point 3 melange
#
h10 = h9
h3 = f*h10 + x7*h2
T3k = PropsSI_nh3_h20('T','P',p3,'H',h3,y3)
T3 = T3k - 273.15
#
#
#
Qevap =   mp*(h7 - h6)
Qcond = mp*(h3 - h4)
Wt = mp*(1-f)*(h1 - h2)
Wp = mp*(h5 - h4)
Wtot = Wt - Wp
eta_c   = Wtot/Qevap
Cp = 4190.0
Cmin1 = mpw*Cp
T3p = T2p - Qevap/Cmin1
h6x = PropsSI_nh3_h20('H','Q',0,'P',p6,y6)
T6xK = PropsSI_nh3_h20('T','Q',0,'P',p6,y6)
T6x = T6xK - 273.15
# latent
#
qb = mp*(h7 - h6x)
T3p = T2p - qb/Cmin1
#
# sensible
#    print('Impossible')
qa = mp*(h6x - h6)
T4p = T3p - qa/Cmin1
if T3p < T6x or T4p < T6:
    print('impossible')
else:
    print('ok')
print('Qevap = ' + '{:5.2f}'.format(Qevap/1000)  + ' kJ/kg')
print('Qcond = ' + '{:5.2f}'.format(Qcond/1000)  + ' kJ/kg')
print('Wtot = ' + '{:5.2f}'.format(Wtot/1000)  + ' kJ/kg')
print('Eout = ' + '{:5.2f}'.format((Qcond+Wtot)/1000)  + ' kJ/kg')
print('eta = '+ '{:5.3f}'.format(eta_c))
print ('\t','T\t\t','p\t\t\t','y\t\t','h')
# point 1
print ('1\t',T1,'\t','%.2f'%(p1/1000),'\t','%.3f'%(y1),'\t','%.2f'%(h1/1000))
# point 2
print ('2\t','%.1f'%T2,'\t','%.2f'%(p2/1000),'\t','%.3f'%(y2),'\t','%.2f'%(h2/1000))
# point 3
print ('3\t','%.1f'%T3,'\t','%.2f'%(p3/1000),'\t','%.3f'%(y3),'\t','%.2f'%(h3/1000))
# point 4
print ('4\t','%.1f'%T4,'\t','%.2f'%(p4/1000),'\t','%.3f'%(y4),'\t','%.2f'%(h4/1000))
# point 5
print ('5\t','%.1f'%T5,'\t','%.2f'%(p5/1000),'\t','%.3f'%(y5),'\t','%.2f'%(h5/1000))
# point 7
print ('7\t','%.1f'%T7,'\t','%.2f'%(p7/1000),'\t','%.3f'%(y7),'\t','%.2f'%(h7/1000))
# point 8
print ('8\t','%.1f'%T8,'\t','%.2f'%(p8/1000),'\t','%.3f'%(y8),'\t','%.2f'%(h8/1000))


