#coding: utf-8
# Exemple 7.8
#
import  numpy as np
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *

pi = np.pi
epsilon = 0
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
DPa = 50   # Pa
cfm = 700
qls = 0.5       # l/s
qms = qls/1000  # m3/s
qmhr = qms*3600 # m3/hr
nbores = 2
nelbows = 5
nvalves = 3
qloopms = qms/nbores  # flow per loop
h_pac = 3   # en mètres
SDR = 11
tube_loop = 1.0
tube_conn = 1.25
hose_diam = 1.0
Hbore = 125.0
bore_sep = 6.0
socket_fusion = False
#
g = 9.8
peak_h = 8.0
COP = 3.5
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
#
# Connecting pipe
#
Cv = Cv_valve('Ball valves',hose_diam)
Cv_st = 20
if socket_fusion:
    Leq_elbow_head = sing_head_loss('Socket 90',tube_conn)
else:   # butt fusion
    Leq_elbow_head = sing_head_loss('Butt 90',tube_conn)
d1,do = sdr_pipe(tube_conn,SDR)
ep = epsilon/d1
A1 = pi*d1**2/4.0
u1 = qms/A1
L1 = 60            # en mètres
Rew = d1*u1/nuw
Re1 = d1*u1/nu
print (' Reynolds number is   = ' + str(Re1) )
ed = ep/d1
fw = Colebrook(Rew,ed)
f1 = Colebrook(Re1,ed)
corr = f1/fw
DPkPa1 = 100*Sg*(qmhr/Cv)**2
DPkPa2 = 100*Sg*(qmhr/Cv_st)**2
h_valve = (nvalves*DPkPa1)*1000/(rho*g) # m de fluide
h_strain = (DPkPa2)*1000/(rho*g) # m de fluide
h_pipe = f1*(L1/d1)*u1**2/(2*g) # Dp en metres
h_minor = nelbows*f1*(Leq_elbow_head/d1)*u1**2/(2*g) # Dp en metres
h_conn = h_strain + h_pac + h_pipe + h_valve + h_minor
print('h elbows  = {:.3f} m'.format(h_minor))
print('h valve  = {:.3f} m'.format(h_valve))
print('h pipe  = {:.3f} m'.format(h_pipe))
print('h strainer  = {:.3f} m'.format(h_strain))
print('h tot conn  = {:.2f} m'.format(h_conn))
#
# loop
#
if socket_fusion:
    Leq_ubend = sing_head_loss('Socket U-bend',tube_loop)
    Leq_elbow  = sing_head_loss('Socket 90',tube_loop)
    Leq_tee_branch  = sing_head_loss('Socket tee-branch',tube_loop)
    Leq_tee_straight  = sing_head_loss('Socket tee-straight',tube_loop)
else:
    Leq_ubend = sing_head_loss('Butt U-bend',tube_loop)
    Leq_elbow  = sing_head_loss('Butt 90',tube_loop)
    Leq_tee_branch  = sing_head_loss('Butt tee-branch',tube_loop)
    Leq_tee_straight  = sing_head_loss('Butt tee-straight',tube_loop)
d2,do = sdr_pipe(tube_loop,SDR)
ed = epsilon/d2
A2 = pi*d2**2/4.0
u2 = qloopms/A2
L2 = 2*Hbore;            # en mètres
Rew = d2*u2/nuw
Re2 = d2*u2/nu
print ('loop Reynolds number is   = ' + str(Re2) )
fw = Colebrook(Rew,ed)
f2 = Colebrook(Re2,ed)

corr = f2/fw
Lm_loop = Leq_ubend +  Leq_elbow + Leq_tee_straight
h_pipe_loop = f2*(L2/d2)*u2**2/(2*g)  # % Dp en metres
h_minor_loop = f2*(Lm_loop/d2)*u2**2/(2*g)  # % Dp en metres
h_loop = h_pipe_loop + h_minor_loop
h_pipe2 = f2*(bore_sep/d2)*u2**2/(2*g)  # % Dp en metres
h_inv = h_loop + h_pipe2
print('h header  = {:.3f} m'.format(h_pipe2))
print('h loop  = {:.3f} m'.format(h_loop))
print('h field  = {:.3f} m'.format(h_inv))

#
#
h_tot = h_inv + h_conn
print('h total  = {:.3f} m'.format(h_tot))
h_pompe = h_tot*1.15         # security factor
print('h total  + safety margin = {:.3f} m'.format(h_pompe))
Dp_tot = rho*g*h_pompe
mp = rho*qms
Wfluid = mp*g*h_pompe
print ('The fluid power is   = {:.1f} W'.format(Wfluid))
Wpompe = mp*g*h_pompe/0.5
Wele = Wpompe/0.85
Qa = ls_cfm(cfm)
Wfan = Qa*DPa/300
HP = Wele/W_hp()
Wpac = peak_h/COP
rap = Wele/peak_h
rap2 = (Wele+Wfan)/(peak_h+ Wfan/1000)
print ('The electric power is   = {:.1f} W'.format(Wele))
print ('The electric power is   = {:.2f} HP'.format(HP))
print ('The fan power is   = {:.2f} W'.format(Wfan))
print ('The ratio W   pumping/kW of heating = {:.2f} W/kW'.format(rap))
print ('The ratio (with fan) W   pumping/kW of heating = {:.2f} W/kW'.format(rap2))
COP_sys = peak_h/(Wpac+Wele/1000)
print ('The system COP (no Fan) is   = {:.2f}'.format(COP_sys))
COP_sys2 = (peak_h+Wfan/1000)/( Wpac+(Wfan+Wele)/1000)
print ('The system COP (Fan) is   = {:.2f}'.format(COP_sys2))


