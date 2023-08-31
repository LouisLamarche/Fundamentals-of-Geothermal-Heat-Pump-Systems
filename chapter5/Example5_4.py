#
# Example 4.4 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md import *
from CoolProp.CoolProp import *
Rp = 0.08
ro = 0.033/2
ri = 0.027/2
kg = 1.0   # grout
ks = 2.5    # soil
rb = 0.075
xc = rb/3.0
xcp = xc/rb
rop = ro/rb


# line source
#
sigma = (kg-ks)/(kg+ks)
Rbc1,Rac1 = Rb_linesource(kg,ks,rb,ro,xc)
Rbc1 = Rbc1 + Rp/2
Rac1 = Rac1 + 2*Rp
xct = xc/rb
rot = ro/rb
Rac2 = 1/(pi*kg)*(np.log(2*xc/ro)+ sigma*np.log((1 + xct**2)/(1 - xct**2))) +2*Rp
print ('Rb line source UT : ' , Rbc1)
print ('Ra line source UT : ' , Rac1)
#
# calcul multi-pole
    #
z = np.array([xc,-xc])  # pipes coordinates
J = 0
Rbe1,Rae1 = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
print ('Rb multipole 10 UT : ',Rbe1)
print ('Ra multipole 10 UT : ',Rae1)


Rb2ut,Ra2ut = Rb_linesource_2UT(kg,ks,rb,ro,xc,'12-34')
Rb2ut = Rb2ut + Rp/4
Ra2ut = Ra2ut + Rp
print ('Rb line source 2UT : ' , Rb2ut)
print ('Ra line source 2UT : ' , Ra2ut)
#
# calcul multi-pole
    #
z = np.array([0 + xc*1j,xc,0-xc*1j,-xc])  # coordonnées des tuyaux
Rb2utb,Ra2utb = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
print ('Rb multipole 10  2UT : ',Rb2utb)
print ('Ra multipole 10  2UT : ',Ra2utb)