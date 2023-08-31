#
# Example 5.2 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md import *

ro = 0.033/2
ri = 0.027/2
kg = 1.0    # grout
ks = 2.5    # soil
rb = 0.075
Rp = 0.08
confv = ['a','b','c']
xcv = [ro+0.00001,rb/3.0,rb-ro-0.00001]
betaov = [20.10,17.44,21.91]
beta1v = [-0.9447,-0.6052,-0.3796]
for i in range(0,3):
    conf = confv[i]
    sconf = 'configuration ' + conf
    xc = xcv[i]
    Rba1 = Rb_Paul(kg,rb,ro,conf) + Rp/2
    betao = betaov[i]
    beta1 = beta1v[i]
    Sg = betao*(rb/ro)**beta1
    Rba2 = 1.0/(Sg*kg) + Rp/2
    print ('Paul-Remuund ' + sconf + ' : ' ,'%.3f' % Rba1,'%.3f' % Rba2)
    #
    # b) Sharqawi
    #
    Rbb1 = Rb_Sharqawi(kg,rb,ro,xc)  + Rp/2
    Rbb2 = 1/(2*pi*kg)*(-1.49*xc/rb+0.656*np.log(rb/ro)+0.436) + Rp/2
    print  ('Sharqawi ' + sconf + ' : ','%.3f' % Rbb1,'%.3f' % Rbb2)
    #
    # line source
    #
    sigma = (kg-ks)/(kg+ks)
    Rbc1,Rac = Rb_linesource(kg,ks,rb,ro,xc)
    Rbc1 = Rbc1 + Rp/2
    xct = xc/rb
    Rbc2 = 1/(4*pi*kg)*(np.log(rb**2/(2*xc*ro))-sigma*np.log(1 - xct**4)) + Rp/2
    print ('line source ' + sconf + ' : ', '%.3f' % Rbc1,'%.3f' % Rbc2)
    #
    # calcul multi-pole
    #
    J = 1
    z = np.array([xc,-xc])  # pipes coordinates
    Rbd1,Rad = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
    #lamda = 1.0
    lamda = (1.0+2*pi*kg*Rp)/ (1.0-2*pi*kg*Rp)
    num =  (1.0-4.0*sigma*xc**4/(rb**4-xc**4))**2
    den =  lamda*(2*xc/ro)**2 + 1.0 + 16.0*sigma*rb**4*xc**4/(rb**4-xc**4)**2
    c =  1/(4*pi*kg)*num/den
    Rbd2 = Rbc2 - c
    print ('multipole 1 ' + sconf + ' : ',Rbd1,Rbd2)
    J = 10
    Rbe1,Rae = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
    print ('multipole 10 ' + sconf + ' : ',Rbe1)
    err1 =abs(Rba1 - Rbe1)/Rbe1*100
    err2 = abs(Rbb1 - Rbe1)/Rbe1*100
    err3 = abs(Rbc1 - Rbe1)/Rbe1*100
    err4 = abs(Rbd1 - Rbe1)/Rbe1*100
#    print('Case '+conf)
    print ('err PAul ' + sconf + ' : ',err1)
    print('err Shar ' + sconf + ' : ',err2)
    print('err line-source ' + sconf + ' : ',err3)
    print('err multi 1 ' + sconf + ' : ',err4)