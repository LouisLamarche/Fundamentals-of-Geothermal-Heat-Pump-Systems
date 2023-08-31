#coding: utf-8
import  numpy as np
from geothermal_md import *
from hydraulic_md import *
from CoolProp.CoolProp import *
from heat_pump_md import *



pi = np.pi
fluide1 = 'water'
fluide2 = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 20 % volume
fluide3 = 'INCOMP::APG-40%'  # ASHRAE propylene glycl 40 % volume
cas = 1
if cas ==1:
    fluide = fluide1
elif cas ==2:
    fluide = fluide2
else:
    fluide = fluide3
patm = 101.325*1000.0
epsilon = 0.0
T_ref = 5.0
T_refk = T_ref+273.15
rhow =PropsSI('D','T',T_refk,'P',patm,'water')
muw = PropsSI('viscosity','T',T_refk,'P',patm,'water')
nuw = muw/rhow
rho =PropsSI('D','T',T_refk,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk,'P',patm,fluide)
exit()
Sg = rho/rhow
nu = mu/rho
g = 9.81
SDR = 11
#
# heat pupms
gpm1 = 9
gpm2 = 9
gpm3 = 8
hc1 = W_BTUhr(50)
hc2 = W_BTUhr(45)
hc3 = W_BTUhr(40)
deb1 = m3s_gpm(gpm1)
deb2 = m3s_gpm(gpm2)
deb3 = m3s_gpm(gpm3)
Dp1 = 1000*kPa_psi(2.0)  # Pascals
Dp2 = 1000*kPa_psi(2.0)
Dp3 = 1000*kPa_psi(1.5)
HP1 = heat_pump(flowrate=deb1,Dp = Dp1,heating_capacity = hc1)
HP2 = heat_pump(flowrate=deb2,Dp = Dp2,heating_capacity = hc2)
HP3 = heat_pump(flowrate=deb3,Dp = Dp3,heating_capacity = hc3)
debit_total = HP1.flowrate + HP2.flowrate + HP3.flowrate
peak_kW = (HP1.heating_capacity + HP2.heating_capacity+ HP3.heating_capacity)/1000
HPv = [HP1,HP2,HP3]
#
#  initial values
# ground
# heat pump
tube = [1.25,1.25,1.25]
#
# hoze kits
#
Cv1 = Cv_valve('Hoze kit',tube[0])
Cv2 = Cv_valve('Hoze kit',tube[1])
Cv3 = Cv_valve('Hoze kit',tube[2])
qmhr1 = HP1.flowrate*3600  # m3s/hr
qmhr2 = HP2.flowrate*3600  # m3s/hr
qmhr3 = HP3.flowrate*3600  # m3s/hr
DPPa1 = Sg*1e5*(qmhr1/Cv1)**2
DPPa2 = Sg*1e5*(qmhr2/Cv2)**2
DPPa3 = Sg*1e5*(qmhr2/Cv3)**2
h_hoze1 = DPPa1/(rho*g) # m de fluide
h_hoze2 = DPPa2/(rho*g) # m de fluide
h_hoze3 = DPPa3/(rho*g) # m de fluide
hhoze = np.array([h_hoze1,h_hoze1,h_hoze3])
#
#
# heaat pumps head calculations
#
tube_pac = 1.5
h_hpv = np.zeros(4)
di,do = sdr_pipe(tube_pac,SDR)
A1 = pi*di**2/4.0
Lpipe = 4
for i in range(0,3):
    q = HPv[i].flowrate
    u = q/A1
    Re = di*u/nu
    ed = epsilon/di
    f = Colebrook(Re,ed)
    Rew = di*u/nuw
    fw = Colebrook(Rew,ed)
    corr = f/fw
    h_pipe = f*((Lpipe)/di)*u**2/(2*g)
    Dp = HPv[i].Dp   # Pascal
    h_hp = corr*Dp/(rho*g)
    h_hpv[i] = h_hp + h_pipe + corr*(hhoze[i])
    print ('Heat pump number #',i+1)
    print(' h heat pump = ',h_hp)
    print(' h piping = ',h_pipe)
    print(' h hoze kit = ',corr*hhoze[i])
    print(' h total = ',h_hpv[i])

#
# building
#
npipes = 6
h_pipe = np.zeros(npipes)
gpmv = np.array([26,17,8,26,9,18])
Lpipe = np.array([60,80,80,40,80,80])
tube = np.array([2,1.5,1.25,2,1.25,1.5])
Le_string = [['Butt tee-straight'],\
            ['Butt tee-straight','Butt 90'],\
            ['Butt 90'],\
            ['Butt tee-straight','Butt 90'],\
            ['Butt tee-straight','Butt 90'],\
            ['Butt tee-straight']]
Leqv = np.zeros(npipes)
for i in range(0,npipes):
    di,do = sdr_pipe(tube[i],SDR)
    A1 = pi*di**2/4.0
    debiti = m3s_gpm(gpmv[i])
    u1 = debiti/A1
    Re = di*u1/nu
    ed = epsilon/di
    f1 = Colebrook(Re,ed)
    Le = 0
    for Le_s in Le_string[i]:
        Leq = sing_head_loss(Le_s,tube[i])
        Le = Le + Leq
    Leqv[i] = Le
    h_pipe[i] = f1*((Lpipe[i]+Le)/di)*u1**2/(2*g) # Dp en metres
    print('pipe number #',i+1,' L equivalent = ','%.2f' % Le,' h total = ',h_pipe[i])
    ff = ft_m(h_pipe[i])/ ft_m(Lpipe[i]+Le)*100
    print('pipe number #',i+1,' friction loss per 100 ft = ','%.3f' % ff)
    ff2 = 9.81*h_pipe[i]/ (Lpipe[i]+Le)
    print('pipe number #',i+1,' friction loss per 100 ft = ','%.3f' % ff2)
hb1 = sum(h_pipe[1:3]) + h_hpv[2]
hb2 = h_pipe[1]  + h_pipe[5] + h_hpv[1]
hb3 = sum(h_pipe[4:6])   + h_hpv[0]
print ('Path second heat pump=',hb3)
print ('Path third heat pump= ',hb2)
print ('Path last heat pump= ',hb1)
h_building = max([hb1,hb2,hb3])
print('h building = ' ,h_building)
h_head = h_pipe[0] + h_pipe[3]
print('h header = ' ,h_head)
exit()
h_ground = 7.5
h_tot = h_building + h_ground + h_head
mp = rho*debit_total
rend = 0.5
Wele = mp*g*h_tot/rend
rap = Wele/peak_kW
print ('The pumping energy is   = ' + str(Wele)  + ' W')
print ('The ratio W   pumping/kW of heating = ' + str(rap) + ' W/kW')




