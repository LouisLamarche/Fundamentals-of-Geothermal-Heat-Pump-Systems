#coding: utf-8
#
# Example 5.1 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md import *
from CoolProp.CoolProp import *
patm = 101.325*1000.0
fluid = 'Water'
ro = 0.021
t = 0.0038
ri = ro - t
di = 2*ri
kg = 1.7    # grout
ks = 2.5    # soil
rb = 0.075
xc = rb/3.0
conf = 'b'
Vp = 0.4                # l/s
Trefk = 30.0+273
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
mp = Vp*1.0e-3*rhof
kp = 0.4       # conductivité du plastique
# pipe resistance
rcond = np.log(ro/ri)/(2*pi*kp)
Re = 4*mp/(pi*di*mu)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hf = (Nud*kf)/(di)
rconv = 1/(pi*di*hf)
Rp = rcond + rconv
Rp2,rcond2,rconv2 = Rp_fct(mp,ro,ri,kp,Trefk,fluid)
print ('Rp',Rp,Rp2)
# a) Paul
Rba1 = Rb_Paul(kg,rb,ro,'b')+ Rp/2.0
betao = 17.44
beta1 = -0.6052
Sg = betao*(rb/ro)**beta1
Rba2 = 1.0/(Sg*kg) + Rp/2.0
print ('Paul-Remund','%.3f' % Rba1,'%.3f' % Rba2)

#
# b) Sharqawi
#
Rbb1 = Rb_Sharqawi(kg,rb,ro,xc)+ Rp/2.0
Rbb2 = 1/(2*pi*kg)*(-1.49*xc/rb+0.656*np.log(rb/ro)+0.436)+ Rp/2.0
print  ('Sharqawi','%.3f' % Rbb1,'%.3f' % Rbb2)
#
# line source
#
sigma = (kg-ks)/(kg+ks)
Rbc1,Rac = Rb_linesource(kg,ks,rb,ro,xc)
Rbc1 = Rbc1 + Rp/2
xct = xc/rb
Rbc2 = 1/(4*pi*kg)*(np.log(rb**2/(2*xc*ro))-sigma*np.log(1 - xct**4))+ Rp/2.0
print ('line source','%.3f' %  Rbc1,'%.3f' % Rbc2)

#
# calcul multi-pole first porder
#
J = 1
z = np.array([xc,-xc])  #
Rbd1,Rad = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
lamda = (1.0+2*pi*kg*Rp)/ (1.0-2*pi*kg*Rp)
num =  (1.0-4.0*sigma*xc**4/(rb**4-xc**4))**2
den =  lamda*(2*xc/ro)**2 + 1.0 + 16.0*sigma*rb**4*xc**4/(rb**4-xc**4)**2
c =  1/(4*pi*kg)*num/den
Rbd2 = Rbc2 - c
print ('multipole 1 ','%.3f' % Rbd1,'%.3f' % Rbd2)
J = 10  # order of multipole
Rbe1,Rae = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
print ('multipole 10 ',Rbe1)
err1 =abs(Rba1 - Rbe1)/Rbe1
err2 = abs(Rbb1 - Rbe1)/Rbe1
err3 = abs(Rbc1 - Rbe1)/Rbe1
err4 = abs(Rbd1 - Rbe1)/Rbe1
print (err1,err2,err3,err4)