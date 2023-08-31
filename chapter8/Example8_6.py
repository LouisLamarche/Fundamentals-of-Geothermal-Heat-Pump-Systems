#
# one-pipe Ex 8.6
#
import numpy as np
from geothermal_md import *
from heat_pump_md import *
from hydraulic_md import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from tabulate import tabulate
npac = 6
VFD = 0
DEB_L = 1400 # load flow rate (gpm)
debit1 = 0.8 # l/s
DEB_S = gpm_ls(debit1)  # source flow rate (gpm)
Tin = 25
Cp = 4.18   # kJ/kg-K
rho = 1000
mp1 = debit1/1000*rho
CCf1 = mp1*Cp
heat_pump = '..\\data\\wa_ts060.xls'
df = pd.read_excel(heat_pump)
tc,sc,po,eer = WA_heat_pump(F_C(Tin),DEB_S,DEB_L,'cooling',df)
#tc,sc,po,eer  = ts060_cooling(F_C(Tin),DEB_S,DEB_L)
q_evap = W_BTUhr(tc)
cop = W_BTUhr(eer)
Wcomp = q_evap/cop
q_cond = q_evap + Wcomp # kW
DeltaTmax = 6
debit_tot = (npac -1)*q_cond/(Cp*DeltaTmax) # kg/s
for i in range(0,5):
    Tout = Tin + q_cond/CCf1
    print('Tout = ',Tout)
    Tin = (mp1*Tout + (debit_tot-mp1)*Tin)/debit_tot
    print('Tin = ',Tin)
debit_tot = max(debit_tot,4.8)
Tin = 25
Tin2 = 25
Energy = np.zeros(npac)
for i in range(0,npac):
    tc,sc,po,eer  = WA_heat_pump(F_C(Tin),DEB_S,DEB_L,'cooling',df)
    q_evap2 = W_BTUhr(tc)
    cop2 = W_BTUhr(eer)
    Wcomp2 = q_evap2/cop2
    q_cond2 = q_evap2 + Wcomp2 # kW
    Tout = Tin + q_cond/CCf1
    Tout2 = Tin2 + q_cond2/CCf1
    Tinn = (Tout*mp1 + Tin*(debit_tot-mp1))/debit_tot
    Tinn2 = (Tout2*mp1 + Tin2*(debit_tot-mp1))/debit_tot
    print('Tout = ',Tout)
    print('Tin = ',Tinn,Tinn2)
    Tin = Tinn
    Tin2 = Tinn2
    Energy[i] = Wcomp2*1000
my_pipes = []

g = 9.81
fac5 = 0.825
fac4 = 0.65
fac3 = 0.49
fac2 = 0.32
nboucles = npac + 1   # = npac
nbouclesi = 0      # nombre de boucles ou le débit est imposé    = Nin + Nout -1
npipes = npac + 13
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
    mboucle = [np.array([1,8,9,10,11,12,13,14,15,16,17,18,19]),np.array([2,-8]),np.array([3,-10]),np.array([4,-12]),np.array([5,-14]),np.array([6,-16]),np.array([7,-18])]
elif npac == 5:
    mboucle = [np.array([1,7,8,9,10,11,12,13,14,15,16,17,18]),np.array([2,-7]),np.array([3,-9]),np.array([4,-11]),np.array([5,-13]),np.array([6,-15])]
elif npac == 4:
    mboucle = [np.array([1,6,7,8,9,10,11,12,13,14,15,16,17]),np.array([2,-6]),np.array([3,-8]),np.array([4,-10]),np.array([5,-12])]
elif npac ==3:
    mboucle = [np.array([1,5,6,7,8,9,10,11,12,13,14,15,16]),np.array([2,-5]),np.array([3,-7]),np.array([4,-9])]
elif npac == 2:
    mboucle = [np.array([1,4,5,6,7,8,9,10,11,12,13,14,15]),np.array([2,-4]),np.array([3,-6])]
# valeurs initiales
debit1 = 0.0008
debit_tot = debit1*npac
debit2 = debit_tot - debit1
qpipe[0] = debit_tot
if npac == 6:
    tube = [3,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3]
elif npac == 5:
    tube = [3,1,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3]
elif npac == 4:
    tube = [3,1,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3]
elif npac == 3:
    tube = [3,1,1,1,3,3,3,3,3,3,3,3,3,3,3,3]
elif npac == 2:
    tube = [3,1,1,3,3,3,3,3,3,3,3,3,3,3,3]
Lpipe[0] = 0.0
for i in range(0,npac):
    ip = i+1
    qpipe[ip] = debit1
    Lpipe[ip] = 0.0
    Kpipe[ip] = 46.5
Lpipe[0] = 325
Lp = 3
for i in range(0,2*npac-1):
    j = i % 2
    debit = j*debit1 + (j-1)*debit2
    ip = npac + i + 1
    qpipe[ip] = debit
    Lpipe[ip] = Lp
qbi[0] =  debit_tot
Lpipe[npipes-1] = 6*Lp

if npac == 6:
    fac = 1
if npac == 5:
    qpipe[npipes-1] =  debit_tot
    Lpipe[npipes-1] = Lp

    fac = fac5
if npac == 4:
    qpipe[npipes-1] =  debit_tot
    qpipe[npipes-2] =  debit_tot
    Lpipe[npipes-1] = Lp
    Lpipe[npipes-2] = Lp
    fac = fac4
if npac == 3:
    qpipe[npipes-1] =  debit_tot
    qpipe[npipes-2] =  debit_tot
    qpipe[npipes-3] =  debit_tot
    Lpipe[npipes-1] = Lp
    Lpipe[npipes-2] = Lp
    Lpipe[npipes-3] = Lp
    fac = fac3

if npac == 2:
    qpipe[npipes-1] =  debit_tot
    qpipe[npipes-2] =  debit_tot
    qpipe[npipes-3] =  debit_tot
    qpipe[npipes-4] =  debit_tot
    Lpipe[npipes-1] = Lp
    Lpipe[npipes-2] = Lp
    Lpipe[npipes-3] = Lp
    Lpipe[npipes-4] = Lp
    fac = fac2
for i in range(1,npac+1):
    qbi[i] = debit1
SDR = 11
for i in range(0,npipes):
    Dpipe[i] =  sdr_pipe(tube[i],SDR)[0]
Pompe =  [np.array([])] * npipes
Pompe[0] = np.array([-1.21547184e+07, -9.03212269e+05,  1.83702735e+03,  6.83208441e+01])

from hpompe_md import pompe_0014,pompe_1935b,pompe_1207
Circ = pompe_0014()
Pompe[0] = pompe_1935b()
for j in range(0,npac):
    Pompe[j+1] = Circ
if VFD:
    for i in range(0,4):
        Pompe[0][i] = Pompe[0][i]*fac**(i-1)
p =Pompe[0]
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i],Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
Trefk = 300
nu = 3.078/3600/1000
my_network = network(nloopi = nbouclesi,mloop = mboucle,my_pipes = my_pipes,nu = nu)
# calul des débits de conduites
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()
dp0 = ghp[0] - gh[0]
rho = 1000
mp = rho*qpipes[0]
mp1 = rho*qpipes[1]
mp2 = rho*qpipes[2]
mp3 = rho*qpipes[3]
mp4 = rho*qpipes[4]
mp5 = rho*qpipes[5]
mp6 = rho*qpipes[6]
W_pompe = mp*ghp[0]
W_pompe1 = mp1*ghp[1]
W_pompe2 = mp2*ghp[2]
W_pompe3 = mp3*ghp[3]
W_pompe4 = mp4*ghp[4]
W_pompe5 = mp5*ghp[5]
W_pompe6 = mp6*ghp[6]
W_circ = W_pompe1 + W_pompe2 +  W_pompe3 + W_pompe4 + W_pompe5 + W_pompe6
deb = qpipes[0]*1000 # debit pompe  l/s
p_rend = np.array([-0.08 , -0.92, 14.50,20.5])
rend_circ = 0.2
P_circ = W_HP(1/8)
rend = np.polyval(p_rend,deb/fac)/100 # efficiency
P_elec = W_pompe/rend + npac*P_circ
#for i in range(0,npac):
#    Energy[i] = Energy[i] +P_elec
print ('Fluid power = ' + str(W_pompe) + ' Watss')
print ('Pump power = ' + str(P_elec) + ' Watss')
h_tot = ghp[0]/g
h_ft = ft_m(h_tot)
Ql = qpipes[0]*1000
gpm = gpm_ls(Ql)
print ('h = ' + str(h_ft) + ' ft')
print ('h = ' + str(h_tot) + ' m')
print ('gpm = ' + str(gpm) + ' gpm')
print ('flow = ' + str(Ql) + ' l/s')
nb = 6
qloop = Ql/nb/1000
di,do = sdr_pipe(1.0,11)
Aloop = pi*di**2/4.0
uloop = qloop/Aloop
Re_loop = di*uloop/nu
print('Re loop = ',Re_loop)
print ('Pump power = ' + '%.2f' % P_elec  + ' W')
head = ['point','Vp(l/s)', 'h(m)']
my_list = []
for i in range(0,npipes):
    vdot = qpipes[i]*1000
    pt = ['%d' % (i+1),'%.2f' % vdot,'%.2f' % (gh[i]/g)]
    my_list.append(pt)
print(tabulate(my_list,headers= head))

