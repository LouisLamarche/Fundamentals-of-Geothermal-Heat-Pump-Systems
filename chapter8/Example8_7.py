import numpy as np
from matplotlib.pyplot import *
from hydraulic_md import *
from conversion_md import *
#
P = 1.5 # HP
Wn = W_HP(P)
eta_pump = 0.8
Wfluid_nom = 573  # from Exemple 8.1
Wfluid_5hp = 400  # from Exemple 8.1
Wbrake_nom = Wfluid_nom/eta_pump
Wbrake_5hp = Wfluid_5hp/eta_pump
xvfd = 0.9      # from Exemple 8.1
xmot_nom  = Wbrake_nom/Wn
xmot_5hp  = Wbrake_5hp/Wn
etam1 = motor_efficiency(P,xmot_nom,case = 'H')
print ('eta mot (6 Heat Pumps) = ',etam1)
eta_v1 = VFD_efficiency(3,1)
print ('eta VFD (6 Heat Pumps) = ',eta_v1)
eta1 = eta_pump*eta_v1*etam1
print ('eta wire-to-water(6 Heat Pumps) = ',eta1)
print ('x ( 5 HP) = ',xmot_5hp)
etam2 = motor_efficiency(P,xmot_5hp,case = 'H')
print ('eta mot (5 Heat Pumps) = ',etam2)
eta_v2 = VFD_efficiency(3,xvfd)
print ('eta VFD (5 Heat Pumps) = ',eta_v2)
eta2 = eta_pump*eta_v2*etam2
print ('eta wire-to-water(5 Heat Pumps) = ',eta2)
