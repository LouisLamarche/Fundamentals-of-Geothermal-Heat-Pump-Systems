import numpy as np
from geothermal_md import *
from hydraulic_md import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from tabulate import tabulate

my_pipes = []
#


npac = 6
VFD = 1
fac = 0.84
g = 9.81
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
mboucle = [np.array([1,2,8,9,10,11,12,18]),np.array([-2,-8,3,13]),np.array([-3,-9,4,14]),np.array([-4,-10,5,15]),np.array([-5,-11,6,16]),np.array([-6,-12,7,17])]
# valeurs initiales
debit1 = 0.0008
debit_nom = debit1*npac
qpipe[0] = debit_nom
tube = [3,1.,1.,1.,1.,1.,1.,1.,1.25,1.5,2,2,2,2,1.5,1.25,1.,3]
Lpipe[0] = 0.0
Lp = 3
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
nu = 3.078/3600/1000
rho = 1000
my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nu)
# calul des débits de conduites
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()
deb = qpipes[0]*1000  # l/s
mp = rho*qpipes[0]
W_pompe = mp*ghp[0]
rend = np.polyval(p_rend,deb/fac)/100 # flow in l/s
P_elec = W_pompe/rend
Dp1 = rho*gh[1]
Dp2 = rho*gh[2]
Dp3 = rho*gh[3]
Dp4 = rho*gh[4]
Dp5 = rho*gh[5]
Dp6 = rho*gh[6]
Dpm = (Dp1+Dp2+Dp3+Dp4+Dp5+Dp6)/6
print ('pressure diff = ' + str(Dpm/1000 ) + ' kPA')
print ('Fluid power = ' + str(W_pompe) + ' Watss')
print ('Pump power = ' + str(P_elec) + ' Watss')
ght1 = (gh[1] + gh[10] + gh[7] + gh[8] + gh[9] + gh[11])/g
print ('Le diff du trajet 1  ' + str(ght1) + ' m')
ght2 = (gh[12] + gh[2] + gh[10] + gh[11] + gh[8] + gh[9])/g
print ('Le diff du trajet 2  ' + str(ght2) + ' m')
ght3 = (gh[12] + gh[13] + gh[3]  +  gh[10] +  gh[11] + gh[9])/g
print ('Le diff du trajet 3  ' + str(ght3) + ' m')
ght4 = (gh[12] + gh[13] + gh[14] + gh[4]  +  gh[10] +  gh[11])/g
print ('Le diff du trajet 4  ' + str(ght4) + ' m')
ght5 = (gh[12] + gh[13] + gh[14] + gh[15] + gh[5] +   gh[11])/g
print ('Le diff du trajet 5  ' + str(ght5) + ' m')
ght6 = (gh[12] + gh[13] + gh[14] + gh[15] +  gh[16] + gh[6])/g
print ('Le diff du trajet 6  ' + str(ght6) + ' m')
deb2 = 1.25*deb
x = np.linspace(0,deb2,10)
p = Pompe[0]/g
y = np.polyval(p,x/1000)
hpi =ghp[0]/g
y2 = hpi*(x/deb)**2
h_tot = ghp[0]/g
h_ft = ft_m(h_tot)
Ql = qpipes[0]*1000
gpm = gpm_ls(Ql)
print ('h = ' + str(h_ft) + ' ft')
print ('h = ' + str(h_tot) + ' m')
print ('gpm = ' + str(gpm) + ' gpm')
print ('flow = ' + str(Ql) + ' l/s')
print ('Pump power = ' + '%.2f' % P_elec  + ' W')
head = ['point','Vp(l/s)', 'h(m)']
my_list = []
for i in range(0,npipes):
    vdot = qpipes[i]*1000
    pt = '%d' % (i+1),'%.2f' % vdot,'%.2f' % (gh[i]/g)
    my_list.append(pt)
print(tabulate(my_list,headers= head))

