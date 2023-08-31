#coding: utf-8
import numpy as np
from scipy.special import *
from hardy_cross_md import *
from CoolProp.CoolProp import *
import pandas as pd
pi = np.pi
patm = 101.325*1000.0
#
nloopi = 0  # number of loops where the flow rate is imposed
nloop = 3
npipes = 6
# vectors
Dpipe = 0.025*np.ones(npipes)
eppipe = 45e-6*np.ones(npipes)
Lepipe = np.zeros(npipes)
Lpipe = np.array([80.0,15.0,10.0,20.0,12.0,25.0])
Kpipe = np.array([2.0,16.0,8.0,0.0,8.0,0.0])
# network topology
mloop = [np.array([1,3,5]),np.array([2,-3,-4]),np.array([4,-5,6])]
# Pump characterics
Pompe = [np.array([-6e7,-1e6,1.9e3,1.8e2]),np.array([-6e7,-1e6,1.9e3,1.8e2]),np.array([]),np.array([]),np.array([]),np.array([])]
my_pipes = []
for i in range(0,npipes):
    my_pipes.append(pipes(Lpipe = Lpipe[i],Dpipe = Dpipe[i],eppipe = eppipe[i], \
    Lepipe = Lepipe[i],Kpipe = Kpipe[i],Pompe = Pompe[i]))
#  initial values
qbi = np.array([0.001,0.0005,0.00075])
Trefk = 300
fluide = 'Water'
rho =PropsSI('D','T',Trefk,'P',patm,fluide)
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluide)
nu = mu/rho
my_network = network(nloopi = nloopi,mloop = mloop,my_pipes = my_pipes,nu = nu)
qb = my_network.hardy_cross(qbi)
## calcul des pertes de charges dans les conduites ( en m2/s2)
qpipes  = my_network.get_qpipes()
gh,ghp  = my_network.get_gh()
if my_network.converged:
    data =  np.zeros((npipes,3))
    for ip in range(0,npipes):
        data[ip,0] = qpipes[ip]*1000  # l/s
        data[ip,1] = gh[ip]
        data[ip,2] = ghp[ip]
    indd  = ['1','2','3','4','5','6']
    coll = ['q(l/s)','(gh)$_{loss}$(J/kg)','(gh)$_{pump}$(J/kg)']
    df = pd.DataFrame(data,columns = coll,index = indd)
    flag_latex = False
    if flag_latex:
        def f1(x):
            return '%.0f' % x
        def f2(x):
            return '%.3f' % x
        def f3(x):
            return '%.2f' % x
        print(df.to_latex(index=True,formatters=[f3,f2,f2]))
    else:
        print(df)
else:
    print('Solution not converged')
