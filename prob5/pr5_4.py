from geothermal_md import *
import  numpy  as np
from conversion_md import *
from CoolProp.CoolProp import *
patm = 101.325*1000.0

# sol
ks = 2.0
Cs = 2.0e6
als = ks/Cs
alhr = als*3600
alj = alhr*24
To = 7.0
load = W_BTUhr(17.9)
cas = 'b'
if cas == 'a':
    CAP = 17.9
    COP = 3.86
    gpm = 5.5
else:
    CAP = 16.7
    COP = 3.69
    gpm = 3.0
qb = W_BTUhr(CAP)*1000
q = qb*(COP-1)/COP
Wcomp = qb/COP
print('W = ',Wcomp/1000,' kW')
Qele =  W_BTUhr(17.9-CAP)
print('W = ',Wcomp/1000 + Qele,' kW')
# fluide
Vp = ls_gpm(gpm)
fluid = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 40 % volume
Trefk = 0.0 + 273
muf = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
Prf = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
Cpf = PropsSI('Cpmass','T',Trefk,'P',patm,fluid)
rhof = PropsSI('D','T',Trefk,'P',patm,fluid)
mpt = Vp*rhof/1000
CCf = mpt*Cpf
#
# puits
#
nb = 1
rb = 0.075
D = 1.0
SDR = 11
di,do = sdr_pipe(D,SDR)
ro = do/2.0
ri = di/2.0
xc = 0.034
if (xc+ro) > rb or xc < ro:
     print('impossible')
     exit()
kp = 0.4
kc = 1
Rcond = np.log(ro/ri)/(2*pi*kp)
rcond = np.log(ro/ri)/(2*pi*kp)
mp1 = mpt/nb
Re = 4*mp1/(pi*di*muf)
print(Re)
Ac = pi*di**2/4
u = mp1/(rhof*Ac)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Prf/8.)/(1.+12.7*np.sqrt(f/8.)*(Prf**(2./3.)-1))
else:
    Nud = 3.6
    print('Careful laminar')
hf = (Nud*kf)/(di)
rconv = 1/(pi*di*hf)
Rp = rcond + rconv
Rp2,rcond2,rconv2 = Rp_fct(mp1,ro,ri,kp,Trefk,fluid)
print ('Rp = ',Rp,Rp2)
Rconv = 1/(pi*di*hf)
Rp = Rcond + Rconv
Rpc,Rpa = Rb_linesource(kc,ks,rb,ro,xc)
Rpb = Rpc + Rp/2
Rpa = Rpa + 2*Rp
print(Rpb,Rpa)
H = 120
L = nb*H
gam = H/(mp1*Cpf*Rpb)
eta = H/(mp1*Cpf*np.sqrt(Rpb*Rpa))
Rbs = Rpb*eta/np.tanh(eta)
qp = q/L
tf = 10
Fo = alhr*tf/rb**2
DTb = -qp*G_function(Fo)/ks
Tb = To + DTb
Tf = Tb - qp*Rbs
Tfo = Tf + q/(2*CCf)
Tfi = Tf - q/(2*CCf)
print(Tfi,Tfo)