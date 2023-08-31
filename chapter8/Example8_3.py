#coding: utf-8
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from tabulate import tabulate

my_pipes = []
#
rho = 1000
nu = 3.078/3600/1000

npac = 2
VFD = 1
g = 9.81
fac6 = 0.5
fac6 = 0.5
fac5 = 0.912
fac4 = 0.83
fac3 = 0.765
fac2 = 0.72


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
p_rend = np.array([-0.08 , -0.92, 14.50,20.5])
# topologie du réseau
if npac == 5:
    mboucle = [np.array([1,2,7,8,9,10,11,16]),np.array([-2,-7,3,12]),np.array([-3,-8,4,13]),np.array([-4,-9,5,14]),np.array([-5,-10,6,15])]
elif npac == 4:
    mboucle = [np.array([1,2,6,7,8,9,10,14]),np.array([-2,-6,3,11]),np.array([-3,-7,4,12]),np.array([-4,-8,5,13])]
elif npac ==3:
    mboucle = [np.array([1,2,5,6,7,8,9,12]),np.array([-2,-5,3,10]),np.array([-3,-6,4,11])]
elif npac == 2:
    mboucle = [np.array([1,2,4,5,6,7,8,10]),np.array([-2,-4,3,9])]


# valeurs initiales
debit1 = 0.0008
debit_nom = debit1*npac
qpipe[0] = debit_nom
if npac == 5:
    tube = [3,1,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,1.25,3]
    if VFD:
        fac = fac5
    else:
        fac =1
    Deltap_setpoint = 61 # kPa

elif npac == 4:
    tube = [3,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,3]
    if VFD:
        fac = fac4
    else:
        fac =1
    Deltap_setpoint = 62 # kPa
elif npac == 3:
    tube = [3,1,1,1,1,1.25,1.5,2,2,2,2,3]
    if VFD:
        fac = fac3
    else:
        fac =1
    Deltap_setpoint = 65 # kPa
elif npac == 2:
    tube = [3,1,1,1,1.25,1.5,2,2,2,3]
    if VFD:
        fac = fac2
    else:
        fac =1
    Deltap_setpoint = 26.0 # kPa
Lp = 3
Lpipe[0] = 0.0
for i in range(0,npac):
    ip = i+1
    qpipe[ip] = debit1
    Lpipe[ip] = 0.0
    Kpipe[ip] = 56
Lpipe[0] = 325
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
ng = 6 - npac
for i in range(0,ng):
    ii = 2*npac + i
    qpipe[ii] =  debit_nom
    Lpipe[ii] = Lp

i1 = npac
for i in range(1,nboucles):
    ip = i1+i
    qbi[i] = qpipe[ip]
SDR = 11
for i in range(0,npipes):
    Dpipe[i] =  sdr_pipe(tube[i],SDR)[0]
Pompe =  [np.array([])] * npipes
Pompe[0] = np.array([-7.80429020e+07, -1.73892164e+05, -9.46230883e+02,  1.36576649e+02])
if VFD:
    for i in range(0,4):
        Pompe[0][i] = Pompe[0][i]*fac**(i-1)
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))

my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nu)
# calul des débits de conduites
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()


dp_setpoint = Deltap_setpoint*1000/rho  # (gh)setpoint J/kg

dp0 = ghp[0] - gh[0]
print('deltaP = ',dp0)

if dp0 > dp_setpoint:
    my_pipesn = []
    nboucles = npac + 1 # = npac
    nbouclesi = 0      # nombre de boucles ou le débit est imposé    = Nin + Nout -1
    npipes = 2*npac + 7
    # initalisation des vecteurs et des cellules
    qbi = np.zeros(nboucles)
    qpipe = np.zeros(npipes)
    Lpipe = np.zeros(npipes)
    eppipe = np.zeros(npipes)
    Lepipe = np.zeros(npipes)
    Kpipe = np.zeros(npipes)
    Dpipe = np.zeros(npipes)


    if npac == 5:
        mboucle = [np.array([1,17]),np.array([2,7,8,9,10,11,16,-17]),np.array([-2,-7,3,12]),np.array([-3,-8,4,13]),np.array([-4,-9,5,14]),np.array([-5,-10,6,15])]

    elif npac == 4:
        mboucle = [np.array([1,15]),np.array([2,6,7,8,9,10,14,-15]),np.array([-2,-6,3,11]),np.array([-3,-7,4,12]),np.array([-4,-8,5,13])]
    elif npac ==3:
        mboucle = [np.array([1,13]),np.array([2,5,6,7,8,9,12,-13]),np.array([-2,-5,3,10]),np.array([-3,-6,4,11])]
    elif npac == 2:
        mboucle = [np.array([1,11]),np.array([2,4,5,6,7,8,10,-11]),np.array([-2,-4,3,9])]


    debit1 = 0.0008
    debit_nom = debit1*npac
    qpipe[0] = debit_nom

    if npac == 5:
        tube = [3,1,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,1.25,3,1]
    elif npac == 4:
        tube = [3,1,1,1,1,1,1.25,1.5,2,2,2,2,1.5,3,1]
    elif npac == 3:
        tube = [3,1,1,1,1,1.25,1.5,2,2,2,2,3,1]
    elif npac == 2:
        tube = [3,1,1,1,1.25,1.5,2,2,2,3,1]
    Lpipe[0] = 0.0
    for i in range(0,npac):
        ip = i+1
        qpipe[ip] = debit1
        Lpipe[ip] = 0.0
        Kpipe[ip] = 56
    Kpipe[0] = 81
    Kpipe[npipes-1] = 800
    i1 = npac+1
    i2 = npipes - 3
    qpipe[i1] = debit1
    qpipe[i2] = debit1
    Lpipe[i1] = Lp
    Lpipe[i2] = Lp
    ng = npac- 2
    for i in range(0,ng):
        ip = i1 + i + 1
        iq = i2 - i  -1
        qpipe[ip] = (i+2)*debit1
        qpipe[iq] = (i+2)*debit1
        Lpipe[ip] = Lp
        Lpipe[iq] = Lp
    qbi[0] =  debit_nom
    qbi[1] =  0.8*debit_nom
    qpipe[ip+1] =  debit_nom
    qpipe[ip+2] =  debit_nom
    Lpipe[ip+1] = Lp
    Lpipe[ip+2] = Lp
    Lpipe[i2+1] = 5*Lp
    i1 = npac + 4
    for i in range(2,nboucles):
        ip = i1+i
        qbi[i] = qpipe[ip]
    SDR = 11
    for i in range(0,npipes):
        Dpipe[i] =  sdr_pipe(tube[i],SDR)[0]
    Pompe =  [np.array([])] * npipes
    Pompe[0] = np.array([-7.80429020e+07, -1.73892164e+05, -9.46230883e+02,  1.36576649e+02])
    if VFD:
        for i in range(0,4):
            Pompe[0][i] = Pompe[0][i]*fac**(i-1)
    for i in range(0,npipes):
        my_pipesn.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
    my_networkn = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipesn,nu = nu)
    qbn = my_networkn.hardy_cross_DPV(qbi,dp_setpoint)
    ## calcul des pertes de charges dans les conduites ( en m2/s2)
    qpipes  = my_networkn.get_qpipes()
    gh,ghp  = my_networkn.get_gh_DPV(dp_setpoint)



deb = qpipes[0]*1000
mp = rho*qpipes[0]
W_pompe = mp*ghp[0]
rend = np.polyval(p_rend,deb/fac)/100 # flow in l/s
P_elec = W_pompe/rend
print ('Fluid power = ' + str(W_pompe) + ' Watss')
print ('Pump power = ' + str(P_elec) + ' Watss')
Dp1 = rho*gh[1]
Dp2 = rho*gh[2]
Dp3 = rho*gh[3]
Dp4 = rho*gh[4]
Dp5 = rho*gh[5]
Dpm = (Dp1+Dp2+Dp3+Dp4+Dp5)/5
print ('pressure diff = ' + str(Dpm/1000 ) + ' kPA')
print(qpipes*1000)
print(ghp-gh)
if npac ==5:
    ght1 = (gh[1] + gh[6] + gh[7] + gh[8] + gh[9]   + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[11] + gh[2] + gh[7] + gh[8] + gh[9] + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght2) + ' m')
    ght3 = (gh[11] + gh[12] + gh[3]  +  gh[8] + gh[9]+ gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght3) + ' m')
    ght4 = (gh[11] + gh[12] + gh[13] + gh[4]  +  gh[9]  + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght4) + ' m')
    ght5 = (gh[11] + gh[12] + gh[13] + gh[14] + gh[5]  + gh[10])/g
    print ('Le diff du trajet 1  ' + str(ght5) + ' m')
elif npac == 4:
    ght1 = (gh[1] + gh[5] + gh[6] + gh[7]  + gh[8] + gh[9])/g
    print ('Le diff du trajet 1  ' + str(ght1) + ' m')
    ght2 = (gh[10] + gh[2] + gh[6] + gh[7] +  gh[8]+ gh[9])/g
    print ('Le diff du trajet 2  ' + str(ght2) + ' m')
    ght3 = (gh[10] + gh[11] + gh[3]  +  gh[7]  +  gh[8]  + gh[9])/g
    print ('Le diff du trajet 3  ' + str(ght3) + ' m')
    ght4 = (gh[10] + gh[11] + gh[12] + gh[4]  +  gh[8]+ gh[9])/g
    print ('Le diff du trajet 4  ' + str(ght4) + ' m')
elif npac ==3 :
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
elif npac ==2 :
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

ecri =0
if ecri:
    if VFD == 0:
        f = open("U:/documents/notes_ang/fig/chapter8/pac52.txt", "w")
    else:
        f = open("U:/documents/notes_ang/fig/chapter8/pac5b2.txt", "w")
    ss = str(deb) + "," + str(hpi)  + ","+  str(ghp[0]) + ","+  str(qp[1]*1000) + ","+  str(gh[1]) +"\n"
    f.write( ss)      # str() converts to string
    ss = str(p[0]) + "," + str(p[1]) + "," +  str(p[2]) + "," + str(p[3]) + "\n"
    f.write( ss)      # str() converts to string
    f.close()
h_tot = ghp[0]/g
h_ft = ft_m(h_tot)
Ql = qpipes[0]*1000
gpm = gpm_ls(Ql)
print ('h = ' + str(h_ft) + ' ft')
print ('h = ' + str(h_tot) + ' m')
print ('gpm = ' + str(gpm) + ' gpm')
print ('flow = ' + str(Ql) + ' l/s')
print ('Pump power = ' + '%.2f' % P_elec  + ' W')
my_list = []
head = ['point','Vp(l/s)', 'h(m)']
for i in range(0,npipes):
    vdot = qpipes[i]*1000
    pt = ['%d' % (i+1),'%.2f' % vdot,'%.2f' % (gh[i]/g)]
    my_list.append(pt)
print(tabulate(my_list,headers= head))
