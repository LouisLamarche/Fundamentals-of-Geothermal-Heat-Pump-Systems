#coding: utf-8
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from tabulate import tabulate

#

npac = 5
VFD = 0
g = 9.81
fac6 = 0.5
fac5 = 0.912
fac4 = 0.83
fac3 = 0.765
fac2 = 0.705
nboucles = npac   # = npac
nbouclesi = 0      # nombre de boucles ou le débit est imposé    = Nin + Nout -1
npipes = 3*nboucles - 1
npipes = 2*npac + 6
# initalisation des vecteurs et des cellules
qbi = np.zeros(nboucles)
qpipe = np.zeros(npipes)
Lpipe = np.zeros(npipes)
eppipe = np.zeros(npipes)
Lepipe = np.zeros(npipes)
Kpipe = np.zeros(npipes)
Dpipe = np.zeros(npipes)
# topologie du réseau
if npac == 6:
    mboucle = [np.array([1,7,13,14,15,16,17,18]),np.array([2,8,-3,-13]),np.array([3,9,-4,-14]),np.array([4,10,-5,-15]),np.array([5,11,-6,-16]),np.array([6,12,-7,-17])]
elif npac == 5:
    mboucle = [np.array([1,6,11,12,13,14,15,16]),np.array([2,7,-3,-12]),np.array([3,8,-4,-13]),np.array([4,9,-5,-14]),np.array([5,10,-6,-15])]
elif npac == 4:
    mboucle = [np.array([1,5,9,10,11,12,13,14]),np.array([2,6,-3,-11]),np.array([3,7,-4,-12]),np.array([4,8,-5,-13])]
elif npac ==3:
    mboucle = [np.array([1,4,7,8,9,10,11,12]),np.array([2,5,-3,-10]),np.array([3,6,-4,-11])]
elif npac == 2:
    mboucle = [np.array([1,3,5,6,7,8,9,10]),np.array([2,4,-3,-9])]
# valeurs initiales
debit1 = 0.0008
debit_nom = debit1*npac
qpipe[0] = debit_nom
if npac == 6:
    tube = [3,1.,1.,1.,1.,1.,1.,1.,1.25,1.5,2,2,2,2,1.5,1.25,1.,3]
elif npac == 5:
    tube = [3,1,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,1.25,3]
elif npac == 4:
    tube = [3,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,3]
elif npac == 3:
    tube = [3,1,1,1,1,1.25,1.5,2,2,2,2,3]
elif npac == 2:
    tube = [3,1,1,1,1.25,1.5,2,2,2,3]
Lpipe[0] = 0.0
for i in range(0,npac):
    ip = i+1
    qpipe[ip] = debit1
    Lpipe[ip] = 0.0
    Kpipe[ip] = 56
Lpipe[0] = 325
Lp = 3
i1 = npac+1
i2 = npipes - 2
i3 = npipes - 1
qpipe[i1] = debit1
qpipe[i2] = debit1
qpipe[i3] = debit_nom
Lpipe[i1] = Lp
Lpipe[i2] = Lp
Lpipe[i3] = Lp*5
ng = npac-2
for i in range(0,ng):
    ip = i1 + i + 1
    iq = i2 - i  -1
    qpipe[ip] = (i+2)*debit1
    qpipe[iq] = (i+2)*debit1
    Lpipe[ip] = Lp
    Lpipe[iq] = Lp
qbi[0] =  debit_nom
if npac == 6:
    if VFD:
        fac = fac6
    else:
        fac =1
if npac == 5:
    qpipe[ip+1] =  debit_nom
    Lpipe[ip+1] = Lp
    if VFD:
        fac = fac5
    else:
        fac =1
if npac == 4:
    qpipe[ip+1] =  debit_nom
    qpipe[ip+2] =  debit_nom
    Lpipe[ip+1] = Lp
    Lpipe[ip+2] = Lp
    if VFD:
        fac = fac4
    else:
        fac =1
if npac == 3:
    qpipe[ip+1] =  debit_nom
    qpipe[ip+2] =  debit_nom
    qpipe[ip+3] =  debit_nom
    Lpipe[ip+1] = Lp
    Lpipe[ip+2] = Lp
    Lpipe[ip+3] = Lp
    if VFD:
        fac = fac3
    else:
        fac =1
if npac == 2:
    qpipe[ip+2] =  debit_nom
    qpipe[ip+3] =  debit_nom
    qpipe[ip+4] =  debit_nom
    qpipe[ip+5] =  debit_nom
    Lpipe[ip+2] = Lp
    Lpipe[ip+3] = Lp
    Lpipe[ip+4] = Lp
    Lpipe[ip+5] = Lp
    if VFD:
        fac = fac2
    else:
        fac =1

i1 = npac
for i in range(1,nboucles):
    ip = i1+i
    qbi[i] = qpipe[ip]
SDR = 11
for i in range(0,npipes):
    Dpipe[i] =  sdr_pipe(tube[i],SDR)[0]
Pompe =  [np.array([])] * npipes
from hpompe_md import pompe_1207
Pompe[0] = pompe_1207()
if VFD:
    for i in range(0,4):
        Pompe[0][i] = Pompe[0][i]*fac**(i-1)
p =Pompe[0]
my_pipes = []
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
Trefk = 300

nu = 3.078/3600/1000
my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nu)
qb = my_network.hardy_cross(qbi)
#
# flow rate calaculation
#
qpipes  = my_network.get_qpipes()
## head losses
gh,ghp  = my_network.get_gh()
dp0 = ghp[0] - gh[0]
rho = 1000
mp = rho*qpipes[0]
W_pompe = mp*ghp[0]
deb = qpipes[0]*1000 # debit pompe  l/s
p_rend = np.array([-0.08 , -0.92, 14.50,20.5])
rend = np.polyval(p_rend,deb/fac)/100 # efficiency
P_elec = W_pompe/rend
print ('Fluid power = ' + str(W_pompe) + ' Watss')
print ('Pump power = ' + str(P_elec) + ' Watss')
flag_aff = True
if npac ==6 and flag_aff:
    ghp1 = gh[1]
    print ('Le diff de la pompe 1 est de = ' + str(ghp1) + ' kPA')
    ghp2 = gh[2]
    ghp3 = gh[3]
    ghp4 = gh[4]
    ghp5 = gh[5]
    ghp6 = gh[6]
    print ('Le diff de la pompe 2 est de = ' + str(ghp2) + ' kPA')
    print ('Le diff de la pompe 3 est de = ' + str(ghp3) + ' kPA')
    print ('Le diff de la pompe 4 est de = ' + str(ghp4) + ' kPA')
    print ('Le diff de la pompe 5 est de = ' + str(ghp5) + ' kPA')
    print ('Le diff de la pompe 6 est de = ' + str(ghp6) + ' kPA')
    ghm = (ghp1+ghp2+ghp3+ghp4+ghp5+ghp6)/6
    print ('Le diff capteur = ' + str(ghm ) + ' kPA')
    ght1 = (ghp1 + gh[10] + gh[7] + gh[8] + gh[9] + gh[11])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[12] + ghp2 + gh[10] + gh[11] + gh[8] + gh[9])/g
    print ('Le diff du trajet 2  ' + str(ght2) + ' m')
    ght3 = (gh[12] + gh[13] + ghp3  +  gh[10] +  gh[11] + gh[9])/g
    print ('Le diff du trajet 3  ' + str(ght3) + ' m')
    ght4 = (gh[12] + gh[13] + gh[14] + ghp4  +  gh[10] +  gh[11])/g
    print ('Le diff du trajet 4  ' + str(ght4) + ' m')
    ght5 = (gh[12] + gh[13] + gh[14] + gh[15] + ghp5 +   gh[11])/g
    print ('Le diff du trajet 5  ' + str(ght5) + ' m')
    ght6 = (gh[12] + gh[13] + gh[14] + gh[15] +  gh[16] + ghp6)/g
    print ('Le diff du trajet 6  ' + str(ght6) + ' m')
elif npac == 5 and flag_aff:
    ghp1 = gh[1]
    print ('Le diff de la pompe 1 est de = ' + str(ghp1) + ' kPA')
    ghp2 = gh[2]
    ghp3 = gh[3]
    ghp4 = gh[4]
    ghp5 = gh[5]
    print ('Le diff de la pompe 2 est de = ' + str(ghp2) + ' kPA')
    print ('Le diff de la pompe 3 est de = ' + str(ghp3) + ' kPA')
    print ('Le diff de la pompe 4 est de = ' + str(ghp4) + ' kPA')
    print ('Le diff de la pompe 5 est de = ' + str(ghp5) + ' kPA')
    ghm = (ghp1+ghp2+ghp3+ghp4+ghp5)/5
    print ('Le diff capteur = ' + str(ghm ) + ' kPA')
    ght1 = (ghp1 + gh[6] + gh[7] + gh[8] + gh[9]   + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[11] + ghp2 + gh[7] + gh[8] + gh[9] + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght2) + 'm')
    ght3 = (gh[11] + gh[12] + ghp3  +  gh[8] + gh[9]+ gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght3) + ' m')
    ght4 = (gh[11] + gh[12] + gh[13] + ghp4  +  gh[9]  + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght4) + ' m')
    ght5 = (gh[11] + gh[12] + gh[13] + gh[14] + ghp5  + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght5) + ' m')
elif npac ==4 and flag_aff:
    ghp1 = gh[1]
    print ('Le diff de la pompe 1 est de = ' + str(ghp1) + ' kPA')
    ghp2 = gh[2]
    ghp3 = gh[3]
    ghp4 = gh[4]
    print ('Le diff de la pompe 2 est de = ' + str(ghp2) + ' kPA')
    print ('Le diff de la pompe 3 est de = ' + str(ghp3) + ' kPA')
    print ('Le diff de la pompe 4 est de = ' + str(ghp4) + ' kPA')
    ghm = (ghp1+ghp2+ghp3+ghp4)/4
    print ('Le diff capteur = ' + str(ghm ) + ' kPA')
#    print ('Le diff capteur = ' + str(ghp4 ) + ' kPA')
    ght1 = (ghp1 + gh[5] + gh[6] + gh[7]  + gh[8] + gh[9])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[10] + ghp2 + gh[6] + gh[7] +  gh[8]+ gh[9])/g
    print ('Le diff du trajet 2  ' + str(ght2) + ' m')
    ght3 = (gh[10] + gh[11] + ghp3  +  gh[7]  +  gh[8]  + gh[9])/g
    print ('Le diff du trajet 3  ' + str(ght3) + ' m')
    ght4 = (gh[10] + gh[11] + gh[12] + ghp4  +  gh[8]+ gh[9])/g
    print ('Le diff du trajet 4  ' + str(ght4) + ' m')
elif npac ==3 and flag_aff:
    ghp1 = gh[1]
    print ('Le diff de la pompe 1 est de = ' + str(ghp1) + ' kPA')
    ghp2 = gh[2]
    ghp3 = gh[3]
    ghm = (ghp1+ghp2+ghp3)/3
    print ('Le diff capteur = ' + str(ghm ) + ' kPA')
    print ('Le diff de la pompe 2 est de = ' + str(ghp2) + ' kPA')
    print ('Le diff de la pompe 3 est de = ' + str(ghp3) + ' kPA')

    ght1 = (ghp1 + gh[4] + gh[5] + gh[6] + gh[7]  + gh[8])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[9] + ghp2 + gh[5] + gh[6] + gh[7]   + gh[8])/g
    print ('Le diff du trajet 1  ' + str(ght2) + ' m')
    ght3 = (gh[9] + gh[10] + ghp3 + gh[6] + gh[7]  + gh[8])/g
    print ('Le diff du trajet 1  ' + str(ght3) + ' m')
elif npac ==2 and flag_aff:
    ghp1 = gh[1]
    print ('Le diff de la pompe 1 est de = ' + str(ghp1) + ' kPA')
    ghp2 = gh[2]
    ghm = (ghp1+ghp2)/2
    print ('Le diff capteur = ' + str(ghm ) + ' kPA')
    print ('Le diff de la pompe 2 est de = ' + str(ghp2) + ' kPA')
    ght1 = (ghp1 + gh[3] + gh[4]+ gh[5]+ gh[6]   + gh[7])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[8] + ghp2 + gh[4] + gh[5] + gh[6] + gh[7])/g
    print ('Le diff du trajet 2  ' + str(ght2) + ' m')
h_tot = ghp[0]/g
h_ft = ft_m(h_tot)
Ql = qpipes[0]*1000
gpm = gpm_ls(Ql)
print ('h = ' + str(h_ft) + ' ft')
print ('h = ' + str(h_tot) + ' m')
print ('gpm = ' + str(gpm) + ' gpm')
print ('flow = ' + str(Ql) + ' l/s')
print ('Fluid power = ' + '%.2f' % W_pompe  + ' W')
print ('Pump power = ' + '%.2f' % P_elec  + ' W')
head = ['point','Vp(l/s)', 'h(m)']
my_list = []
for i in range(0,npipes):
    vdot = qpipes[i]*1000
    pt = ['%d' % (i+1),'%.2f' % vdot,'%.2f' % (gh[i]/g)]
    my_list.append(pt)
print(tabulate(my_list,headers= head))

X = ght1/h_tot
y = fct_lemire(5/6,X)
q = y*415
