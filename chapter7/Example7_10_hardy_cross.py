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
nbouclesi = 1       # nomber of pipes where flow is imposed
nboucles = 4
npipes = 12
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
mboucle = [np.array([2,6,-5,-9]),np.array([3,7,-6,-10]), \
        np.array([4,8,-7,-11]),np.array([1,5,9,10,11,12])]
#  initial values
debit1 = 0.00025
socket_fusion = False
tube_conn =1.5
tube_loop = 1.0
qtot = debit1*4
if socket_fusion:
    Leq_elbow = sing_head_loss('Socket 90',tube_conn)
else:   # butt fusion
    Leq_elbow = sing_head_loss('Butt 90',tube_conn)        #
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
L_head = 30
L_loop = 250
L_seg = 6.0
Dpipe[0] = D15
Dpipe[1] = D125
Dpipe[2] = D125
Dpipe[3] = D1
Dpipe[4] = D1
Dpipe[5] = D1
Dpipe[6] = D1
Dpipe[7] = D1
Dpipe[8] = D1
Dpipe[9] = D125
Dpipe[10] = D125
Dpipe[11] = D15
qpipe[0] = qtot
qpipe[1] = 3*debit1
qpipe[2] = 2*debit1
qpipe[3] = debit1
qpipe[4] = debit1
qpipe[5] = debit1
qpipe[6] = debit1
qpipe[7] = debit1
qpipe[8] = 3*debit1
qpipe[9] = 2*debit1
qpipe[10] = debit1
qpipe[11] = qtot
Lpipe[0] = L_head
Lpipe[1] = L_seg
Lpipe[2] = L_seg
Lpipe[3] = L_seg
Lpipe[4] = L_loop
Lpipe[5] = L_loop
Lpipe[6] = L_loop
Lpipe[7] = L_loop
Lpipe[8] = L_seg
Lpipe[9] = L_seg
Lpipe[10] = L_seg
Lpipe[11] = L_head
nelbows = 5
Lepipe[0] =   nelbows*Leq_elbow
Lepipe[5] = Leq_ubend   + Leq_tee_straight + Leq_elbow
Lepipe[6] = Leq_ubend   + Leq_tee_straight + Leq_elbow
Lepipe[7] = Leq_ubend   + Leq_tee_straight + Leq_elbow
Pompe = [np.array([]),np.array([]),np.array([]),np.array([]),np.array([]), \
    np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])]
my_pipes = []
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
qbi[0] = qpipe[1]
qbi[1] = qpipe[2]
qbi[2] = qpipe[3]
qbi[3] = qpipe[0]
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