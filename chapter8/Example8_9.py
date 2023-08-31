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
Sg = rho/rhow
nu = mu/rho
g = 9.81
SDR = 11
#
# heat pupms
ton1 = 6
ton2 = 5
gpm1 = 16
gpm2 = 13
hc1 = W_BTUhr(ton1*12000)
hc2 = W_BTUhr(ton2*12000)
deb1 = m3s_gpm(gpm1)
deb2 = m3s_gpm(gpm2)
Dp1 = 1000*kPa_psi(4.0)  # Pascals
Dp2 = 1000*kPa_psi(3.4)
HP1 = heat_pump(flowrate=deb1,Dp = Dp1,heating_capacity = hc1)
HP2 = heat_pump(flowrate=deb2,Dp = Dp2,heating_capacity = hc2)
debit_total = 2*HP1.flowrate + 2*HP2.flowrate
peak_kW = (2*HP1.heating_capacity + 2*HP2.heating_capacity)/1000
HPv = [HP1,HP1,HP2,HP2]

#
#  initial values
# ground


# heat pump
tube = [1.25,1.25]
#
# hoze kits
#
Cv1 = Cv_valve('Hoze kit',tube[0])
Cv2 = Cv_valve('Hoze kit',tube[1])
qmhr1 = HP1.flowrate*3600  # m3s/hr
qmhr2 = HP2.flowrate*3600  # m3s/hr
DPPa1 = Sg*1e5*(qmhr1/Cv1)**2
DPPa2 = Sg*1e5*(qmhr2/Cv2)**2
h_hoze1 = DPPa1/(rho*g) # m de fluide
h_hoze2 = DPPa2/(rho*g) # m de fluide
#
# strainers
#
Cv_str1 = Cv_valve('Y-strainer',tube[0])
Cv_str2 = Cv_valve('Y-strainer',tube[1])
DPPas1 = Sg*1e5*(qmhr1/Cv_str1)**2
DPPas2 = Sg*1e5*(qmhr2/Cv_str2)**2
h_str1 = DPPas1/(rho*g) # m de fluide
h_str2 = DPPas2/(rho*g) # m de fluide
hhoze = np.array([h_hoze1,h_hoze1,h_hoze2,h_hoze2])
hstr = np.array([h_str1,h_str1,h_str2,h_str2])
#
#  ball valves
#
Cv_val1 = Cv_valve('Ball valves',tube[0])
Cv_val2 = Cv_valve('Ball valves',tube[1])
DPPas1 = Sg*1e5*(qmhr1/Cv_val1)**2
DPPas2 = Sg*1e5*(qmhr2/Cv_val2)**2
h_val1 = DPPas1/(rho*g) # m de fluide
h_val2 = DPPas2/(rho*g) # m de fluide
hvalve = np.array([h_val1,h_val1,h_val2,h_val2])


#
# heaat pumps head calculations
#
tube_pac = 1.5
h_hpv = np.zeros(4)
di,do = sdr_pipe(tube_pac,SDR)
A1 = pi*di**2/4.0
Lpipe = 4
for i in range(0,4):
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
    h_hpv[i] = h_hp + h_pipe + corr*(hhoze[i] + hstr[i] + 2*hvalve[i])
    print ('Heat pump number #',i+1)
    print(' h heat pump = ',h_hp)
    print(' h piping = ',h_pipe)
    print(' h hoze kit = ',corr*hhoze[i])
    print(' h strainer = ',corr*hstr[i])
    print(' h valve = ',2*corr*hvalve[i])
    print(' h total = ',h_hpv[i])

#
# building
#
h_pipe = np.zeros(8)
gpmv = np.array([58,42,26,13,58,16,32,45])
Lpipe = np.array([12,10,66,50,26,10,66,50])
tube = np.array([3,3,2,1.5,3,1.5,2,3])
Le_string = [['Butt tee-straight','Butt 90'],['Butt tee-straight'],\
            ['Butt tee-straight','Butt 90','Butt 90','Butt reducer'],['Butt 90','Butt 90'],[],\
            ['Butt tee-straight','Butt 90'],['Butt tee-straight','Butt 90'],\
            ['Butt tee-straight','Butt 90']]
Leqv = np.zeros(8)
for i in range(0,8):
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
hb1 = sum(h_pipe[1:4]) + h_hpv[3]
hb2 = sum(h_pipe[1:3])  + h_pipe[7] + h_hpv[2]
hb3 = sum(h_pipe[1:2])  + sum(h_pipe[6:8]) + h_hpv[1]
hb4 =  sum(h_pipe[5:8]) + h_hpv[0]
print ('Path first heat pump= ',hb4)
print ('Path second heat pump=',hb3)
print ('Path third heat pump= ',hb2)
print ('Path last heat pump= ',hb1)
h_building = max([hb1,hb2,hb3,hb4])
print('h building = ' ,h_building)
h_head = h_pipe[0] + h_pipe[4]
print('h header = ' ,h_head)

h_ground = 7.5
h_tot = h_building + h_ground + h_head
mp = rho*debit_total
rend = 0.5
Wele = mp*g*h_tot/rend
rap = Wele/peak_kW
print ('The pumping energy is   = ' + str(Wele)  + ' W')
print ('The ratio W   pumping/kW of heating = ' + str(rap) + ' W/kW')




