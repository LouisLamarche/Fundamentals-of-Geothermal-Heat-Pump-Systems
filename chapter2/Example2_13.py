import numpy as np
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
# exemple 2.5
patm = 101.325*1000.0
Ti = 22.0
To = 30.0
Vol = 72.0
phii = 0.3
phio = 0.8
DT = Ti-To
Tm = (Ti+To)/2.0
Tmk = Tm+273
Cp = PropsSI('Cpmass','T',Tmk,'P',patm,'air')
rho = PropsSI('D','T',Tmk,'P',patm,'air')
mp = 0.5*Vol/3600*rho
K1 = mp*Cp
q = K1*DT
print ('%.2f' % q + ' W')
Tik = Ti+273.15
Tok = To+273.15
psin = PropsSI('P','T',Tik,'Q',0,'water')
psout = PropsSI('P','T',Tok,'Q',0,'water')
pwin = phii*psin
pwout = phio*psout
wi1 = 0.62198*pwin/(patm-pwin)
wo1 = 0.62198*pwout/(patm-pwout)
wi = HAPropsSI('W','T',Tik,'R',phii,'P',patm)
wo = HAPropsSI('W','T',Tok,'R',phio,'P',patm)
hf = PropsSI('H','T',Tmk,'Q',0,'water')
hg = PropsSI('H','T',Tmk,'Q',1,'water')
hfg = hg - hf
q1 = mp*(wo-wi)*hfg
print (' Latent loads are:')
# exemple de nombre formaté à 2 décimales
print ('%.2f' % q1 + ' W')
q2 = mp*(wo1-wi1)*hfg
print (' Latent loads are:')
# exemple de nombre formaté à 2 décimales
print ('%.2f' % q2 + ' W')