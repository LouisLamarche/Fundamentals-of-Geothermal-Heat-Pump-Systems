#coding: utf-8
#   fonctions utiles en hydraulique
#       Colebrook(Re,eD,flag_turb=0)
import numpy as np
from geothermal_md import *


pente = 12000
pi = np.pi
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])              # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])   # nombre d'heures ?coul?es apr?s chaque mois
facteur = 5e4
class borefield:

    def __init__(self,params = [],ground = [],borehole = []):
        self.params  = params
        self.ground = ground
        self.Tp =  0
        self.eta =  0
        self.Rbs =  0
        self.borehole  = borehole

    def Compute_L_ashrae(self):
        qa = self.params.qa
        qh_ch = self.params.qh_heating
        qh_cl = self.params.qh_cooling
        qm_ch = self.params.qm_heating
        qm_cl = self.params.qm_cooling
        flag_Tp = self.params.flag_Tp
        flag_inter = self.params.flag_inter
        flag_Ch = self.params.flag_Ch
        als = self.ground.alsoil
        ks = self.ground.ksoil
        To = self.ground.To
        Rb = self.borehole.Rb
        rb = self.borehole.rb
        CCf = self.borehole.CCf
        d = self.borehole.dist
        nannees = self.params.n_years
        nbloc = self.params.n_bloc
        nrings = self.params.n_rings
        nx = self.borehole.nx
        ny = self.borehole.ny
        nb = nx*ny
        CCu = CCf/nb
        Tfo_ch = self.params.Tfo_heating
        Tfo_cl = self.params.Tfo_cooling
        Tfi_ch = Tfo_ch - qh_ch/(CCf)
        Tfi_cl = Tfo_cl - qh_cl/(CCf)
        Tf_ch = (Tfi_ch+Tfo_ch)/2.0
        Tf_cl = (Tfi_cl+Tfo_cl)/2.0
        DeltaTch  = To - Tf_ch
        DeltaTcl = To - Tf_cl
        alhr = als*3600
        alj = alhr*24
        t1 = nannees*8760
        t2 = t1 + 730 # nombre d'annees + 1 mois
        t3 = t2 + nbloc  # nombre d'annees + 1 mois + bloc horaire
        Fo1 = alhr*(t3-t1)/rb**2
        Fo2 = alhr*(t3-t2)/rb**2
        Fof = alhr*t3/rb**2
        Gf = G_function(Fof)
        G1 = G_function(Fo1)
        G2 = G_function(Fo2)
        #
        # calcul des r?sistances de sol
        #
        Rann = (Gf-G1)/ks
        Rmen = (G1-G2)/ks
        Rhor = G2/ks
        #
        alj = alhr*24
        #
        Fsc = 1.0
        Tp = 0
        ok = False
        Li = 500.0
        delta = 0.1
        compt = 1
        qa_cl = min(qa,0)
        qa_ch = max(qa,0)

        while not ok :
            Hi = Li/nb
            if flag_inter:
                Ra = self.borehole.Rint
                eta = Hi/(CCu*np.sqrt(Rb*Ra))
                Rbs = Rb*eta/np.tanh(eta)
                self.eta = eta
                self.Rbs = Rbs
            else:
                Rbs = Rb
                self.Rbs = Rbs
            Tpch = Tp*qa_ch/(qa_ch+0.0001)
            Tpcl = Tp*qa_cl/(qa_cl+0.0001)
            Hch1 = (qa_ch*Rann)/(DeltaTch - Tpch)
            Hch2 = (qm_ch*Rmen)/(DeltaTch - Tpch)
            Hch4 = (qh_ch*Rbs)/(DeltaTch - Tpch)
            Hch3 = (qh_ch*Rhor*Fsc)/(DeltaTch- Tpch)
            H_ch = Hch1+Hch2+Hch3+Hch4
            # climatisation

            Hcl1 = qa_cl*Rann/(DeltaTcl-Tpcl)
            Hcl2 = qm_cl*Rmen/(DeltaTcl-Tpcl)
            Hcl3 = qh_cl*Rhor*Fsc/(DeltaTcl-Tpcl)
            Hcl4 = qh_cl*Rbs/(DeltaTcl-Tpcl)
            H_cl = Hcl1+Hcl2+Hcl3+Hcl4
            L = max(H_ch,H_cl)
            H = L/nb
            qap1 = qa/L  # W/m
            if flag_Tp == 'ASHRAE':
                if flag_Ch:
                    Tpn = Tp_ashrae(nx,ny,d,qap1,alj,ks,nannees,nrings,H = H)
                else:
                    Tpn = Tp_ashrae(nx,ny,d,qap1,alj,ks,nannees,nrings)
            elif flag_Tp == 'ILS':
                Tpn = Tp_ils(nx,ny,d,qap1,alj,ks,nannees)
            elif flag_Tp == 'FLS':
                Tpn = Tp_fls(nx,ny,d,qap1,alj,ks,nannees,H)
            else:
                Tpn = Tp_Bernier(nx,ny,d,qap1,alj,ks,nannees,H)
            if (abs(L-Li)) < delta:
                ok = True
            else:
                compt = compt + 1
                if compt > 1000:
                    print ('divergence')
                    ok = True
                Tp = Tpn
                Li = L
        self.Tp = Tp
        self.L_he = H_ch
        self.L_co = H_cl
        self.Tfi_he = Tfi_ch
        self.Tfi_co = Tfi_cl
        return L




    def Compute_L_eed(self):


        als = self.ground.alsoil
        alhr = als*3600
        alj = alhr*24
        ks = self.ground.ksoil
        To = self.ground.To
        Rb = self.borehole.Rb
        rb = self.borehole.rb
        CCf = self.borehole.CCf
        d = self.borehole.dist
        nannees = self.params.n_years
        zo = self.params.zo
        nbloc = self.params.n_bloc
        nx = self.borehole.nx
        ny = self.borehole.ny
        nb = nx*ny
        CCu = CCf/nb
        Tfo_ch = self.params.Tfo_heating
        Tfo_cl = self.params.Tfo_cooling
        mois_ini = self.params.init_month
        q_sol_mois = self.params.q_months
        mois_ch = np.argmax(q_sol_mois)
        mois_cl = np.argmin(q_sol_mois)
        qa = np.mean(q_sol_mois)
        if qa > 0:
            nannees_ch = nannees
            nannees_cl = 0
            if mois_cl - mois_ini < 0:
                nannees_cl = 1
        else:
            nannees_cl = nannees
            nannees_ch = 0
            if mois_ch - mois_ini < 0:
                nannees_ch = 1
        qh_ch = self.params.qh_heating
        qh_cl = self.params.qh_cooling
        flag_inter = self.params.flag_inter
        #
        #
        Tfi_ch = Tfo_ch - qh_ch/(CCf)
        Tfi_cl = Tfo_cl - qh_cl/(CCf)
        Tf_ch = (Tfi_ch+Tfo_ch)/2.0
        Tf_cl = (Tfi_cl+Tfo_cl)/2.0
        Tf_ch = (Tfi_ch+Tfo_ch)/2.0
        Tf_cl = (Tfi_cl+Tfo_cl)/2.0
        DeltaTch  = To - Tf_ch
        DeltaTcl = To - Tf_cl
        #
        nb = nx*ny
        CCu = CCf/nb
        Hi = 100.0
        #
        # d?but du calcul
        #
        ntot_ch = int(nannees_ch*12 + mois_ch + 1 - mois_ini)  + 1
        ntot_cl = int(nannees_cl*12 + mois_cl + 1 - mois_ini)  + 1
        tf_ch = (ntot_ch-1)*730.0 + nbloc  # temps final en heures
        tf_cl = (ntot_cl-1)*730.0 + nbloc  # temps final en heures
        tf = max(tf_ch,tf_cl)

        ok = False
        delta = 0.5
        compt = 1
        compt_max = 20
        while not ok :
            rr = rb/Hi
            zob = zo/Hi
            ts = Hi**2/(9*alhr)         # temps caract?ristique Eskilson en heures
            dx = d/Hi
            dy = d/Hi
            if flag_inter:
                Ra = self.borehole.Rint
                eta = Hi/(CCu*np.sqrt(Rb*Ra))
                Rbs = Rb*eta/np.tanh(eta)
                self.eta = eta
                self.Rbs = Rbs
            else:
                Rbs = Rb
                self.Rbs = Rbs
            #
            # g?n?ration du champ de capteur rectangulaire
            #
            nb = nx*ny
            zt = np.zeros([nb,2])
            k = 0
            for i1  in range(0,nx):
                x = i1*dx
                for i2 in range(0,ny):
                    y = i2*dx
                    zt[k] = [x,y]
                    k = k+1
            #
            ng = 80
            t1 = 1.0/ts         # temps minimum 1 heure
            t2 = (1.5*tf)/ts    # temps maximum 1.5 fois tfinal
            ttt = np.logspace(np.log10(t1),np.log10(t2),ng)      # On calcule ng valeurs de la fonction g qu'on interpole par la suite
            ntt = len(ttt)
            nt = ntt+1
            tt = np.zeros(nt)
            g1 = np.zeros(nt)
            tt[1:nt] = ttt[0:ntt]
            for i in range (0,nt):
                g1[i] = compute_g_function(zt,tt[i],rbb = rr,zob = zob)
            #
            q1 = q_sol_mois[mois_ini]
            Foo = tf_ch/ts
            summ  = q1*np.interp(Foo,tt,g1)
            ti = 0
            for ia in range(1,ntot_ch-1):
               nm = (ia + mois_ini) % 12   # On trouve le mois dans l'annee
               ti = ti + 730.0
               Foo = (tf_ch-ti)/ts
               summ = summ + (q_sol_mois[nm]-q1)*np.interp(Foo,tt,g1)
               q1 = q_sol_mois[nm]
            Foo = nbloc/ts
            # ajout du pulse horaire
            simm = summ + (qh_ch-q1)*np.interp(Foo,tt,g1)
            somm = simm/(2*pi*ks)+qh_ch*Rbs
            Hch_Es = somm/DeltaTch
            # clim
            q1 = q_sol_mois[mois_ini]
            Foo = tf_cl/ts
            summ  = q1*np.interp(Foo,tt,g1)
            ti = 0
            for ia in range(1,ntot_cl-1):
               nm = (ia + mois_ini) % 12   # On trouve le mois dans l'annee
               ti = ti + 730.0
               Foo = (tf_cl-ti)/ts
               summ = summ + (q_sol_mois[nm]-q1)*np.interp(Foo,tt,g1)
               q1 = q_sol_mois[nm]
            Foo = nbloc/ts
            # ajout du pulse horaire
            summ = summ + (qh_cl-q1)*np.interp(Foo,tt,g1)
            somm = summ/(2*pi*ks)+qh_cl*Rbs
            Hcl_Es = somm/DeltaTcl
            H_Es = max(Hch_Es,Hcl_Es)
            Hin = H_Es/nb
            if abs(Hin-Hi) < delta:
                ok = True
            else:
                Hi = Hin
                compt = compt + 1
                if compt> compt_max:
                    err =1
                    ok = True
        self.Rbs = Rbs
        self.L_he = Hch_Es
        self.L_co = Hcl_Es
        self.Tfi_he = Tfi_ch
        self.Tfi_co = Tfi_cl

        return H_Es


class ashrae_params:

    """

    Attributes
    ----------
    Lpipe : float
        pipe length (in meters).
    Dpipe : float
        pipe diameter (in meters).

    """


    def __init__(self, qa=0, qm_heating=0, qm_cooling=0, qh_heating=0, qh_cooling=0,Tfo_heating = 0,\
                Tfo_cooling = 25,n_years = 20, n_bloc = 4,n_rings = 3,flag_Tp = 'ASHRAE',flag_inter = False,flag_Ch = False):
        self.qa = float(qa)
        self.qm_heating = float(qm_heating)
        self.qm_cooling = float(qm_cooling)
        self.qh_heating = float(qh_heating)
        self.qh_cooling = float(qh_cooling)
        self.Tfo_heating = float(Tfo_heating)
        self.Tfo_cooling = float(Tfo_cooling)
        self.n_years = int(n_years)
        self.n_rings = int(n_rings)
        self.n_bloc = int(n_bloc)
        self.flag_Tp = flag_Tp
        self.flag_inter = flag_inter
        self.flag_Ch = flag_Ch


class eed_params:

    """

    Attributes
    ----------
    Lpipe : float
        pipe length (in meters).
    Dpipe : float
        pipe diameter (in meters).

    """


    def __init__(self, init_month = 0, q_months =[], qh_heating=0, qh_cooling=0, \
        Tfo_heating = 0,Tfo_cooling = 25,n_years = 20, n_bloc = 4,zo = 4,flag_inter = False):
        self.init_month = int(init_month)
        self.q_months = q_months
        self.qh_heating = float(qh_heating)
        self.qh_cooling = float(qh_cooling)
        self.Tfo_heating = float(Tfo_heating)
        self.Tfo_cooling = float(Tfo_cooling)
        self.zo = float(zo)
        self.n_years = int(n_years)
        self.n_bloc = int(n_bloc)
        self.flag_inter = flag_inter


class borehole:

    """

    Attributes
    ----------

    """

    def __init__(self, nx=1, ny = 1, dist = 6.0,rb = 0.07,Rb = 0.1, CCf=2000,Rint = 0.4):
        self.nx = int(nx)
        self.ny = int(ny)
        self.CCf = float(CCf)
        self.Rb = float(Rb)
        self.Rint = float(Rint)
        self.rb = rb
        self.dist = dist

class ground:

    """

    Attributes
    ----------
    ksoil: float
        ground conductivity (W/m-K).
    alsoil: float
        ground diffusivity (m2/s).
    To : float
        ground temperature (C ).

    """
    def __init__(self, ksoil=1, alsoil = 1e-6, To = 10):
        self.ksoil = float(ksoil)
        self.alsoil = float(alsoil)
        self.To = float(To)
