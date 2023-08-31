#coding: utf-8
#
# 20 sept 2018
#
import numpy as  np
from CoolProp.CoolProp import *
fluide = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 20 % volume
patm = 101.325*1000.0
Trefk = 4 + 273.15
rho =PropsSI('D','T',Trefk,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',Trefk,'P',patm,fluide)

#
# data
#
patm = 101.325*1000.0
#
debit = 0.001     # debit volumique
mpch = debit*rho
mpcl = mpch
Cch = mpch*Cp
Ccl = mpcl*Cp
Qch = 30000.0
Qcl = 30000.0
#
# (i)
#
def COPcooling(T):
    return 8 - 0.18*T + 0.001*T**2
def COPheating(T):
    return 3.5 + 0.06*T - 0.0005*T**2
Tinch = 4.0                     # Tin first PAC
COPch  = COPheating(Tinch)
Qinch = (COPch-1)/COPch*Qch     # evaporator
Pcha = Qch/COPch                 # 1st HP power
Toutcha = Tinch -Qinch/Cch      # temperature à la sortie
Tincl = 4.0                # Tin 2eme  PAC
COPcl  = COPcooling(Tincl)
Qincl = (COPcl+1)/COPcl*Qcl    # puissance au condenseur  2eme PAC
Pcla = Qcl/COPcl                 # consommation de la 2eme  PAC
Toutcla = Tincl + Qincl/Ccl       # temperature à la sortie
Ta = (mpch*Toutcha + mpcl*Toutcla)/(mpch+mpcl) # temperature moyenne
Pa = Pcha+Pcla            # consommation totale a)
#
# 2 ème  configuration
#
Tinch = 4.0
COPch  = COPheating(Tinch)
Qinch = (COPch-1)/COPch*Qch
Pchb = Qch/COPch
Toutchb = Tinch -Qinch/Cch
Tincl = Toutchb              # Tin 2eme  PAC
COPcl  = COPcooling(Tincl)
Qincl = (COPcl+1)/COPcl*Qcl
Pclb = Qcl/COPcl
Toutclb = Tincl + Qincl/Ccl         # temperature à la sortie
Tib = Tincl
Tb = Toutclb
Pb = Pchb+Pclb
#
# 3 ème  configuration
#
Tincl = 4.0
COPcl  = COPcooling(Tincl)
Qincl = (COPcl+1)/COPcl*Qcl
Pclc = Qcl/COPcl
Toutclc = Tincl + Qincl/Ccl
Tinch = Toutclc
COPch  = COPheating(Tinch)
Qinch = (COPch-1)/COPch*Qch
Pchc = Qch/COPch
Toutchc = Tinch - Qinch/Cch
Tc = Toutchc
Pc = Pchc+Pclc
print ('reponse 2a')
print ('%.2f' % (Pa/1000), 'kW')
print ('%.2f' % (Pb/1000), 'kW')
print ('%.2f' % (Pc/1000), 'kW')
print ('reponse 2b')
print ('%.2f' % (Ta), ' C')
print ('%.2f' % (Tb), ' C')
print ('%.2f' % (Tc), ' C')
