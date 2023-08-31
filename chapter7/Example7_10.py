#coding: utf-8
# Exemple 7.10
#
import  numpy as np
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *

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
qls = 1.0       # l/s
q_nom = 1.0
qms = qls/1000  # m3/s
qmhr = qms*3600 # m3/hr
nbores = 4
nelbows = 5
nvalves = 3
qloopms = qms/nbores  # flow per loop
h_pac = 3*(qls/q_nom)**2   # en mètres
SDR = 11
tube_loop = 1.0
tube_conn = 1.5
hose_diam = 1.0
Hbore = 125.0
bore_sep = 6.0#
socket_fusion = False
#
ed = 0.0
g = 9.8
peak_h = 16.0
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
#
# Connecting pipe
#
Cv = Cv_valve('Ball valves',hose_diam)
Cv_st = 20  # m3/hr/bar
if socket_fusion:
    Leq_elbow = sing_head_loss('Socket 90',tube_conn)
else:   # butt fusion
    Leq_elbow = sing_head_loss('Butt 90',tube_conn)
d1,do = sdr_pipe(tube_conn,SDR)
A1 = pi*d1**2/4.0
u1 = qms/A1
L1 = 60            # en mètres
Rew = d1*u1/nuw
Re1 = d1*u1/nu
print ('Reynolds number is   = ' + str(Re1) )
fw = Colebrook(Rew,ed)
f1 = Colebrook(Re1,ed)
corr = f1/fw

DPkPa1 = 100*Sg*(0.5*qmhr/Cv)**2
DPkPa2 = 100*Sg*(0.5*qmhr/Cv_st)**2
h_valve = (nvalves*DPkPa1 + DPkPa2)*1000/(rho*g) # m de fluide
#h_valve = nvalves*DPkPa*1000/(rho*g) # m de fluide
h_pipe = f1*(L1/d1)*u1**2/(2*g) # Dp en metres
h_minor = nelbows*f1*(Leq_elbow/d1)*u1**2/(2*g) # Dp en metres
h_conn = h_pac + h_pipe + h_valve + h_minor
h_hp = h_pac + h_valve
print('h elbows  = {:.3f} m'.format(h_minor))
print('h tot conn  = {:.2f} m'.format(h_conn))

#
# loop
#
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
d_loop,do = sdr_pipe(tube_loop,SDR)
A_loop = pi*d_loop**2/4.0
u_loop = qloopms/A_loop
L_loop = 2*Hbore;            # en mètres
Rew = d_loop*u_loop/nuw
Re_loop = d_loop*u_loop/nu
print ('loop Reynolds number is   = ' + str(Re_loop) )
f_loop = Colebrook(Re_loop,ed)
Lm_loop = Leq_ubend +  Leq_elbow + Leq_tee_straight
h_pipe_loop = f_loop*(L_loop/d_loop)*u_loop**2/(2*g)  # % Dp en metres
h_minor_loop = f_loop*(Lm_loop/d_loop)*u_loop**2/(2*g)  # % Dp en metres
h_loop = h_pipe_loop + h_minor_loop
print('h loop  = {:.3f} m'.format(h_loop))
# pipe 2 , pipe 11
d2,do = sdr_pipe(1.25,SDR)
q2 = 3*qloopms
A2 = pi*d2**2/4.0
u2 = q2/A2
Re2 = d2*u2/nu
print ('loop Reynolds number (pipe 2-11) is   = ' + str(Re2) )
f2 = Colebrook(Re2,ed)
h_2 = f2*(bore_sep/d2)*u2**2/(2*g)  # % Dp en metres
# pipe 3 , pipe 10
d3,do = sdr_pipe(1.25,SDR)
q3 = 2*qloopms
A3 = pi*d3**2/4.0
u3 = q3/A3
Re3 = d3*u3/nu
print ('loop Reynolds number (pipe 3-10) is   = ' + str(Re3) )
f3 = Colebrook(Re3,ed)
h_3 = f3*(bore_sep/d3)*u3**2/(2*g)  # % Dp en metres
# pipe 4 , pipe 9
d4,do = sdr_pipe(tube_loop,SDR)
q4 = qloopms
A4 = pi*d4**2/4.0
u4 = q4/A4
Re4 = d4*u4/nu
print ('loop Reynolds number (pipe 4-9) is   = ' + str(Re4) )
f4 = Colebrook(Re4,ed)
h_4 = f4*(bore_sep/d4)*u4**2/(2*g)  # % Dp en metres
print('h pipe 2  = {:.3f} m'.format(h_2))
print('h pipe 3  = {:.3f} m'.format(h_3))
print('h pipe 4  = {:.3f} m'.format(h_4))

hin1 = h_2 + h_3 + h_4
print('h  path 2-3-4-8  = {:.3f} m'.format(hin1))
hin2 = h_2 + h_3 + h_2
print('h  path 2-3-7-11  = {:.3f} m'.format(hin2))

hinv = 0
tube = [1.0,1.25,1.25]
for i in range(0,nbores-1):
    i1 = i+1.0
    x = i1/nbores
    qms2 = qms*x          #  débit qui passe dans la section i
    D,DD =  sdr_pipe(tube[i],SDR)
    A = pi*D**2/4.0
    u = qms2/A
    Re = D*u/nu
    f = Colebrook(Re,0.0)
    h = f*(bore_sep/D)*u**2/(2*g) #  en metres
    hinv = hinv + 2*x*h
#    print(hinv)
heq1 = hin1+h_loop
print('h  field path 2-3-4-8  = {:.2f} m'.format(heq1))
heq2 = hin2+h_loop
print('h  field path 2-3-7-11  = {:.2f} m'.format(heq2))
heq3 = hinv+h_loop
print('h  field Eq 7.26  = {:.2f} m'.format(heq3))
#
#
h_tot1 = heq1 + h_conn
h_tot2 = heq2 + h_conn
h_tot3 = heq3 + h_conn
print('h  total path 2-3-4-8  = {:.2f} m'.format(h_tot1))
print('h  total path 2-3-7-11  = {:.2f} m'.format(h_tot2))
print('h  total Eq 7.26  = {:.2f} m'.format(h_tot3))
h_pompe = h_tot3*1.15         # security factor
#
#
#
#
print('h total  + safety margin (Eq 7.26) = {:.3f} m'.format(h_pompe))

Dp_tot = rho*g*h_pompe
mp = rho*qms
Wfluid = mp*g*h_pompe
print ('The fluid power is   = {:.1f} W'.format(Wfluid))
Wpompe = mp*g*h_pompe/0.5
Wele = Wpompe/0.85
HP = Wele/W_hp()
Wpac = peak_h/COP
rap = Wele/peak_h
print ('The electric power is   = {:.1f} W'.format(Wele))
print ('The electric power is   = {:.2f} HP'.format(HP))
print ('The ratio W   pumping/kW of heating = {:.2f} W/kW'.format(rap))
COP_sys = peak_h/(Wpac+Wele/1000)
print ('The system COP  is   = {:.2f}'.format(COP_sys))


