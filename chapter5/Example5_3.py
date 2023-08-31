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
xcv = [ro+0.0005,rb/3.0,rb-ro-0.0005]
confv = ['a','b','c']

# line source
#
sigma = (kg-ks)/(kg+ks)
for i in range(0,3):
    xc = xcv[i]
    conf = confv[i]
    sconf = 'configuration ' + conf

    Rbc1,Rac1 = Rb_linesource(kg,ks,rb,ro,xc)
    Rac1 = Rac1 + 2*Rp
    xct = xc/rb
    rot = ro/rb
    Rac2 = 1/(pi*kg)*(np.log(2*xc/ro)+ sigma*np.log((1 + xct**2)/(1 - xct**2))) +2*Rp
    print ('line source' + sconf + ' : ' , '%.3f' % Rac1,'%.3f' % Rac2)
    #
    # calcul multi-pole
    #
    J = 1
    z = np.array([xc,-xc])  # pipes coordinates
    Rbd1,Rad1 = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
    lamda = (1.0+2*pi*kg*Rp)/ (1.0-2*pi*kg*Rp)
    num = 1/(pi*kg)*((ro/(2*xc))**2*(1+ 4*sigma*xct**2/(1 - xct**4))**2)
    den = lamda - (ro/(2*xc))**2 + 2.0*sigma*rot**2*(1+xct**4)/(1-xct**4)**2
    Rad2 = Rac1 - num/den
    print ('multipole 1 ' + sconf + ' : ' ,'%.3f' % Rad1,'%.3f' % Rad2)
    J = 10
    Rbe1,Rae1 = Rb_multipole(kg,ks,rb,ro,Rp,J,z)
    Rab1 = np.arccosh((2*xc**2 - ro**2)/(ro**2))/(2*pi*kg) + 2*Rp
    print ('shape factor  ' + sconf + ' : ' ,'%.3f' % Rab1)

    print ('multipole 10 ' + sconf + ' : ' ,'%.3f' % Rae1)
    err1 = abs(Rac1 - Rae1)/Rae1 *100
    print('err line-source ' + sconf + ' : ',err1)
    err2 = abs(Rad1 - Rae1)/Rae1*100
    print('err  multi-pole order 1  ' + sconf + ' : ',err1)
    err3 = abs(Rab1 - Rae1)/Rae1*100
    print('err shape factor ' + sconf + ' : ',err3)
