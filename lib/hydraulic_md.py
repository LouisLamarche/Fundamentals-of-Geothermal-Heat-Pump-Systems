#coding: utf-8
#   fonctions utiles en hydraulique
#       Colebrook(Re,eD,flag_turb=0)
import numpy as np
from scipy.optimize import newton
from scipy.integrate import *
from scipy.special import *
from collections import *
from conversion_md import *
import pandas as pd




pente = 31000
pi = np.pi
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])              # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])   # nombre d'heures ?coul?es apr?s chaque mois
facteur = 5e4
class network:

    def __init__(self,nloopi=0,mloop=[],my_pipes=[],nu =  1e-6):
        self.nloopi  = int(nloopi)
        self.mloop  = mloop
        self.my_pipes  = my_pipes
        self.converged = False
        self.nu = nu  # fluid viscosity m2/s
    def  hardy_cross(self,qbi):
        mboucle = self.mloop
        nboucles = len(mboucle)
        nbouclesi = nboucles - self.nloopi
        pipe = self.my_pipes
        npipes = len(pipe)

        def  gct_HC(qb):
            qp = np.zeros(npipes)
            y = np.zeros(nbouclesi)
            den = np.zeros(nbouclesi)
            for ib in range(0,nboucles):
                pipesv = mboucle[ib]
                npi = len(pipesv)
                for ip in range(0,npi):
                    jup = pipesv[ip]
                    sg = np.sign(jup)
                    jp = abs(jup)-1
                    qp[jp] = qp[jp] + sg*qb[ib]
            for ib in range(0,nbouclesi):
                pipesv = mboucle[ib]
                npi = len(pipesv)
                somgh = 0
                somd = 0
                for ip in range(0,npi):
                    jup = pipesv[ip]
                    sg = np.sign(jup)
                    jp = np.abs(jup)-1
                    L = pipe[jp].Lpipe
                    D = pipe[jp].Dpipe
                    ep = pipe[jp].eppipe
                    q = qp[jp]
                    sq = np.sign(q)
                    Re = 4*np.abs(q)/(pi*D*self.nu)
                    epD = ep/D
                    ppompe = pipe[jp].Pompe
                    ll = len(ppompe)
                    if ll:
                        dppompe = np.array([ppompe[0]*3,ppompe[1]*2,ppompe[2]])
                    else:
                        dppompe = np.array([])
                    fn = Colebrook(Re,epD)
                    gh = (fn*(L/D+ pipe[jp].Lepipe/D)+pipe[jp].Kpipe)*8*q**2/(pi**2*D**4)
                    gh_pompe = np.polyval(ppompe,q)
                    somgh = somgh + sq*sg*gh - sq*sg*gh_pompe
                    d_pomp = np.polyval(dppompe,q)
                    if q ==0:
                        somd == 1
                    else:
                        somd = somd+np.abs(2*gh/q - d_pomp)
                y[ib] = somgh
                den[ib] = somd
            return y,den
        tol = 1e-6
        count_max = 300
        ok = False
        count = 1
        err = 0
        qb = np.copy(qbi)
        dq = np.zeros(nbouclesi)
        while not ok:
            y,den  = gct_HC(qb)
            ymax = max(abs(y))
            if ymax < tol:
                ok = True
                self.converged = True
            else:
                for ib in range(0,nbouclesi):
                    dq[ib] = -y[ib]/den[ib]
                    qb[ib] = qb[ib] + dq[ib]
                count = count + 1;
                if count > count_max:
                    err = 1
                    ok = True
        qp = np.zeros(npipes)
        for ib in range(0,nboucles):
            pipesv = mboucle[ib]
            npi = len(pipesv)
            for ip in range(0,npi):
                jup = pipesv[ip]
                sg = np.sign(jup)
                jp = abs(jup)-1
                qp[jp] = qp[jp] + sg*qb[ib]
        for jp  in range(0,npipes):
            self.my_pipes[jp].qpipe = qp[jp]
        return qb

    def  hardy_cross_DPV(self,qbi,set_point):
        mboucle = self.mloop
        nboucles = len(mboucle)
        nbouclesi = nboucles - self.nloopi
        pipe = self.my_pipes
        npipes = len(pipe)

        def  gct_HC(qb):
            qp = np.zeros(npipes)
            y = np.zeros(nbouclesi)
            den = np.zeros(nbouclesi)
            for ib in range(0,nboucles):
                pipesv = mboucle[ib]
                npi = len(pipesv)
                for ip in range(0,npi):
                    jup = pipesv[ip]
                    sg = np.sign(jup)
                    jp = abs(jup)-1
                    qp[jp] = qp[jp] + sg*qb[ib]
            for ib in range(0,nbouclesi):
                pipesv = mboucle[ib]
                npi = len(pipesv)
                somgh = 0
                somd = 0
                for ip in range(0,npi):
                    jup = pipesv[ip]
                    sg = np.sign(jup)
                    jp = np.abs(jup)-1
                    L = pipe[jp].Lpipe
                    D = pipe[jp].Dpipe
                    ep = pipe[jp].eppipe
                    q = qp[jp]
                    sq = np.sign(q)
                    Re = 4*np.abs(q)/(pi*D*self.nu)
                    epD = ep/D
                    ppompe = pipe[jp].Pompe
                    ll = len(ppompe)
                    if ll:
                        dppompe = np.array([ppompe[0]*3,ppompe[1]*2,ppompe[2]])
                    else:
                        dppompe = np.array([])
                    fn = Colebrook(Re,epD)
                    if jp == npipes - 1:
                        gh =  (set_point + pente*np.abs(q))
                    else:
                        gh = (fn*(L/D+ pipe[jp].Lepipe/D)+ pipe[jp].Kpipe)*8*q**2/(pi**2*D**4)
                    gh_pompe = np.polyval(ppompe,q)
                    somgh = somgh + sq*sg*gh - sq*sg*gh_pompe
                    d_pomp = np.polyval(dppompe,q)
                    somd = somd+np.abs(2*gh/q - d_pomp)
                y[ib] = somgh
                den[ib] = somd
            return y,den
        tol = 1e-6
        count_max = 200
        ok = False
        count = 1
        err = 0
        qb = np.copy(qbi)
        dq = np.zeros(nbouclesi)
        while not ok:
            y,den  = gct_HC(qb)
            ymax = max(abs(y))
            if ymax < tol:
                ok = True
                self.converged = True
            else:
                for ib in range(0,nbouclesi):
                    dq[ib] = -y[ib]/den[ib]
                    qb[ib] = qb[ib] + dq[ib]
                count = count + 1;
                if count > count_max:
                    err = 1
                    ok = True
        qp = np.zeros(npipes)
        for ib in range(0,nboucles):
            pipesv = mboucle[ib]
            npi = len(pipesv)
            for ip in range(0,npi):
                jup = pipesv[ip]
                sg = np.sign(jup)
                jp = abs(jup)-1
                qp[jp] = qp[jp] + sg*qb[ib]
        for jp  in range(0,npipes):
            self.my_pipes[jp].qpipe = qp[jp]
        return qb


    def  get_qpipes(self):
        pipe = self.my_pipes
        npipes = len(pipe)
        qp = np.zeros(npipes)
        for i in range(0,npipes):
            qp[i] = pipe[i].qpipe
        return qp


    def  get_gh(self):
        pipe = self.my_pipes
        npipes = len(pipe)
        gh= np.zeros(npipes)
        ghp = np.zeros(npipes)
        gh_pompe = np.zeros(npipes)
        for jp in range(0,npipes):
            L = pipe[jp].Lpipe
            D = pipe[jp].Dpipe
            ep = pipe[jp].eppipe
            q = pipe[jp].qpipe
            sq = np.sign(q)
            Re = 4*np.abs(q)/(pi*D*self.nu)
            epD = ep/D
            ppompe = pipe[jp].Pompe
            ll = len(ppompe)
            if ll:
                dppompe = np.array([ppompe[0]*3,ppompe[1]*2,ppompe[3]])
            else:
                dppompe = np.array([])
            fn = Colebrook(Re,epD)
            ghp[jp] = (fn*(L/D+ pipe[jp].Lepipe/D)+pipe[jp].Kpipe)*8*q**2/(pi**2*D**4)
            gh_pompe[jp] = np.polyval(ppompe,q)
            gh[jp] = ghp[jp] - gh_pompe[jp]
        return ghp,gh_pompe
    def  get_gh_DPV(self,set_point):
        pipe = self.my_pipes
        npipes = len(pipe)
        gh= np.zeros(npipes)
        ghp = np.zeros(npipes)
        gh_pompe = np.zeros(npipes)
        for jp in range(0,npipes):
            L = pipe[jp].Lpipe
            D = pipe[jp].Dpipe
            ep = pipe[jp].eppipe
            q = pipe[jp].qpipe
            sq = np.sign(q)
            Re = 4*np.abs(q)/(pi*D*self.nu)
            epD = ep/D
            ppompe = pipe[jp].Pompe
            ll = len(ppompe)
            if ll:
                dppompe = np.array([ppompe[0]*3,ppompe[1]*2,ppompe[3]])
            else:
                dppompe = np.array([])
            fn = Colebrook(Re,epD)
            if jp == npipes - 1:
                ghp[jp] =  (set_point + pente*np.abs(q))
            else:
                ghp[jp] = (fn*(L/D+ pipe[jp].Lepipe/D)+pipe[jp].Kpipe)*8*q**2/(pi**2*D**4)
            gh_pompe[jp] = np.polyval(ppompe,q)
            gh[jp] = ghp[jp] - gh_pompe[jp]
        return ghp,gh_pompe

class pipes:

    """

    Attributes
    ----------
    Lpipe : float
        pipe length (in meters).
    Dpipe : float
        pipe diameter (in meters).

    """

    def __init__(self, Lpipe=0, Dpipe=0, eppipe=0, Lepipe=0, Kpipe=0,qpipe = 0,Pompe=[]):
        self.Lpipe = float(Lpipe)      # Borehole length
        self.Dpipe = float(Dpipe)      # Borehole buried depth
        self.eppipe = float(eppipe)  # Borehole radius
        self.Lepipe = float(Lepipe)      # Borehole x coordinate position
        self.Kpipe = float(Kpipe)      # Borehole y coordinate position
        self.qpipe = float(qpipe)      # Borehole y coordinate position
        self.Pompe = Pompe

def sing_head_loss(type,siz):
    size = '%.2f' % float(siz)
#    size = '%.2f' % siz
    df = pd.read_csv("../data/SI_equivalent_lengths.csv")
    df.fillna(-99,inplace = True)
    df  = df.set_index('type')
    df = df.rename(lambda x: x.strip())
    if type in df.index:
        if size in df.columns:
            val = df.at[type,size]
            if val == -99:
                print('No value given for this size')
        else:
            val = -99
            print('No value given for this size')
    else:
        print('You should choose between these types:')
        print(df.index.values)
        val = -99
    return float(val)
def Cv_valve(type,siz):
    size = '%.2f' % float(siz)
    df = pd.read_csv("../data/SI_Cv_valves.csv")
    df.fillna(-99,inplace = True)
    df  = df.set_index('type')
    df = df.rename(lambda x: x.strip())
    if not isinstance(size, str):
        size = str(float(size))
    if type in df.index:
        if size in df.columns:
            val = df.at[type,size]
            if val == -99:
                print('No value given for this size')
        else:
            val = -99
            print('No value given for this size')
    else:
        print('You should choose between these types:')
        print(df.index.values)
        val = -99
    return float(val)

def Colebrook(Re,eD=0,flag_turb=0):

    def  Cole(f):
        y = 1/np.sqrt(f)+2*np.log10(eD/3.7+2.51/(Re*np.sqrt(f)))
        return y


    # si flag_turb = 1, on force la turbulence m?me si Re < 2300
    if ( Re > 2300 or flag_turb == 1):
        if eD ==0:
            # conduite lisse
             f = (0.790*np.log(Re)-1.64)**(-2)
        else:
            fi = 0.3164/Re**0.25
            f = newton(Cole,fi)
    else:
        if Re <= 0:
            f = 0
        else:
        # laminaire
            f = 64.0/Re
    return f

def f_smooth(Re):
    f = (0.790*np.log(Re)-1.64)**(-2)
    return f


def pump_power(q,h,rho = 1000,eta = 1):  # Compute HP from h ( ft) and q (gpm)
    g = 9.81
    hm = m_ft(h)
    qm = m3s_gpm(q)
    W = qm*rho*g*hm/eta
    HP = HP_W(W)
    return HP




def VFD_efficiency(HP,x_vfd):
    VFD = np.array([[3,0.978856,0.034247,-0.007862],\
        [5,0.977485,0.028413,-0.002733] ,\
        [10,0.978715,0.022227,0.001941],\
        [20,0.984973,0.017545,-0.000475],\
        [30,0.987405,0.015536,-0.005937],\
        [50,0.987910,0.018376,-0.001692],\
        [60,0.971904,0.014537,0.011849],\
        [75,0.991874,0.017897,-0.001301],\
        [100,0.982384,0.012598,0.001405]])
    av = VFD[:,1]
    bv = VFD[:,2]
    cv = VFD[:,3]
    xv = VFD[:,0]
    a2 = np.interp(HP,xv,av)
    b2 = np.interp(HP,xv,bv)
    c2 = np.interp(HP,xv,cv)
    eta =  a2*x_vfd/(b2+x_vfd) + c2*x_vfd
    return eta

def motor_efficiency(HP,x_mot,case = 'H'):
    mot = np.array([[1,1.092165,0.08206,-0.0072],\
        [5,1.223684,0.08467,-0.135186],\
        [10,1.146258,0.045766,-0.110367],\
        [25,1.137209,0.050236,-0.08915],\
        [50,1.088803,0.029753,-0.064058],\
        [75,1.07714,0.029005,-0.04935],\
        [100,1.035294,0.012948,-0.024708]])
    am = mot[:,1]
    bm = mot[:,2]
    cm = mot[:,3]
    xm = mot[:,0]

    mot_max = np.array([[0.196205,3.653654,0.839926],\
        [0.292280,3.368739,0.762471], \
        [0.395895,3.065240,0.674321]])
    def eta_max(cas,P):
        Y = np.log(P)
        if cas == 'H':
            amm = mot_max[0,0]
            bmm = mot_max[0,1]
            cmm = mot_max[0,2]
        elif cas =='M':
            amm = mot_max[1,0]
            bmm = mot_max[1,1]
            cmm = mot_max[1,2]
        elif cas =='L':
            amm = mot_max[2,0]
            bmm = mot_max[2,1]
            cmm = mot_max[2,2]
        eta =  amm*Y/(bmm+Y) + cmm
        return eta
    a2 = np.interp(HP,xm,am)
    b2 = np.interp(HP,xm,bm)
    c2 = np.interp(HP,xm,cm)
    etamax = eta_max(case,HP)
    etam =  a2*x_mot/(b2+x_mot) + c2*x_mot
    eta = etam*etamax
    return eta


def fct_lemire(q,x):
    y = 0.22464 - 0.50812*q + 0.58266*q**2 + 0.7912*q**3 - 0.23834*x \
    + 0.06207*x**2 + 0.15208*x**3 + 2.05530*q*x - 0.15778*q*x**2 \
    - 1.6652*q**2*x -0.17377*q**2*x**2
    return y