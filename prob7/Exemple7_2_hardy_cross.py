#coding: latin-1
from numpy import *
from geothermal_md import *
from scipy.special import *
from collections import namedtuple
from hydraulic_md import *
from CoolProp.CoolProp import *

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

T_refk = 300
rhow =PropsSI('D','T',T_refk,'P',patm,'water')
muw = PropsSI('viscosity','T',T_refk,'P',patm,'water')
nuw = muw/rhow
rho =PropsSI('D','T',T_refk,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk,'P',patm,fluide)
Sg = rho/rhow
nu = mu/rho
g = 9.81
pipes = namedtuple("pipes", "Lpipe Dpipe eppipe Lepipe Kpipe Pompe")
my_pipes = []
#
nbouclesi = 2
nbouclese = 1      # nomber of pipes where flow is imposed
nboucles = nbouclesi + nbouclese
npipes = 8
SDR = 11
D2 = 0.02
D3 = 0.03
qbi = zeros(nboucles)
qpipe = zeros(npipes)
Dpipe = zeros(npipes)
Lpipe = zeros(npipes)
eppipe = zeros(npipes)
Lepipe = zeros(npipes)
Kpipe = zeros(npipes)
# network
mboucle = [array([2,5,-7,-4]),array([3,6,-8,-5]),array([1,4,7,8])]
#  initial values
qtot = 1.2/1000
debit1 = qtot/3
# loop
L_loop =45.974641753487006
L_seg = 8.0
Dpipe[0] = D3
Dpipe[1] = D3
Dpipe[2] = D2
Dpipe[3] = D2
Dpipe[4] = D2
Dpipe[5] = D2
Dpipe[6] = D2
Dpipe[7] = D3
qpipe[0] = qtot
qpipe[1] = 2*debit1
qpipe[2] = debit1
qpipe[3] = debit1
qpipe[4] = debit1
qpipe[5] = debit1
qpipe[6] = debit1
qpipe[7] = 2*debit1
Lpipe[0] = 0
Lpipe[1] = L_seg
Lpipe[2] = L_seg
Lpipe[3] = L_loop
Lpipe[4] = L_loop
Lpipe[5] = L_loop
Lpipe[6] = L_seg
Lpipe[7] = L_seg
Kpipe[4] = 4
Pompe = [array([]),array([]),array([]),array([]),array([]),array([]),array([]),array([])]
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
qbi[0] = qpipe[1]
qbi[1] = qpipe[2]
qbi[2] = qpipe[0]
qb = hardy_cross(qbi,nboucles,nbouclesi,mboucle,my_pipes,npipes,nu)
# calul des débits de conduites
qp = zeros(npipes)
for ib in range(0,nboucles):
    pipesv = mboucle[ib]
    np = len(pipesv)
    for ip in range(0,np):
        jup = pipesv[ip]
        sg = sign(jup)
        jp = abs(jup)-1
        qp[jp] = qp[jp] + sg*qb[ib]
#  head losses
gh,ghp = fct_gh(qp,npipes,my_pipes,nu)
W = qp*gh
hloop1 = (gh[1] + gh[2] + gh[5])/g;print(hloop1)
hloop2 = (gh[1] + gh[4] + gh[7] )/g;print(hloop2)
hloop3 = (gh[3] + gh[6] + gh[7] )/g;print(hloop3)
qp1 = qp[3]*1000;print('Flow in borehole 1=',qp1)
qp2 = qp[4]*1000;print('Flow in borehole 2=',qp2)
qp3= qp[5]*1000;print('Flow in borehole 3=',qp3)
