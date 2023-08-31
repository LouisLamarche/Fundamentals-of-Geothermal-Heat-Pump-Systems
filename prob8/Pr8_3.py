##################################coding: latin-1
import  numpy as np
from geothermal_md import *
from hydraulic_md import *
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
#
nbouclesi = 3       # nomber of  loops where flow is imposed
nboucles = 3
npipes = 7
SDR = 11
D1,do = sdr_pipe(1,SDR)
D125,do = sdr_pipe(1.25,SDR)
D15,do = sdr_pipe(1.5,SDR)
qbi = np.zeros(nboucles)
qpipe = np.zeros(npipes)
Dpipe = np.zeros(npipes)
Lpipe = np.zeros(npipes)
eppipe = np.zeros(npipes)
Lepipe = np.zeros(npipes)
Kpipe = np.zeros(npipes)
# network
mboucle = [np.array([1,2,3,7]),np.array([4,-5,-2]), \
        np.array([5,6,-3])]
#  initial values
gpm1 = 8
gpm2 = 12
gpm3 = 15
debit1 = m3s_gpm(gpm1)
debit2 = m3s_gpm(gpm2) + debit1
debit3 = m3s_gpm(gpm3) + debit2
Dpipe[0] = D15
Dpipe[1] = D125
Dpipe[2] = D125
Dpipe[3] = D1
Dpipe[4] = D1
Dpipe[5] = D1
Dpipe[6] = D1
qpipe[0] = debit3
qpipe[1] = debit3-debit1
qpipe[2] = debit3
qpipe[3] = debit1
qpipe[4] = m3s_gpm(gpm3)
qpipe[5] = debit2
qpipe[6] = debit3
Lpipe[0] = 100
Lpipe[1] = 10
Lpipe[2] = 10
Lpipe[3] = 10
Lpipe[4] = 10
Lpipe[5] = 10
Lpipe[6] = 100
Kpipe[2] = 2
Kpipe[3] = 2
Kpipe[4] = 2
Pompe = [np.array([]),np.array([]),np.array([]),np.array([]),np.array([]), \
    np.array([]),np.array([])]
my_pipes = []
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
qbi[0] = qpipe[0]
qbi[1] = qpipe[2]
qbi[2] = qpipe[5]
my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nu)
# calul des débits de conduites
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()
W = qpipes*gh
h_head = (gh[0] + gh[11])/g
heq = sum(W)/qtot/g - h_head
hloop1 = (gh[1] + gh[2] + gh[3] + gh[7])/g;print('h field path 1=' ,hloop1)
hloop2 = (gh[1] + gh[2] + gh[6] + gh[10])/g;print('h field path 2=' ,hloop2)
hloop3 = (gh[1] + gh[5] + gh[9] + gh[10])/g;print('h field path 3=' ,hloop3)
hloop4 = (gh[4] + gh[8] + gh[9] + gh[10])/g;print('h field path 4=' ,hloop4)
qp1 = qpipes[4]*1000;print('Flow in borehole 1=',qp1)
qp2 = qpipes[5]*1000;print('Flow in borehole 2=',qp2)
qp3= qpipes[6]*1000;print('Flow in borehole 3=',qp3)
qp4 = qpipes[7]*1000;print('Flow in borehole 4=',qp4)