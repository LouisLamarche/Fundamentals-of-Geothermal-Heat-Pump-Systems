#coding: utf-8
import numpy as np
import pandas as pd
from geothermal_md import *
from conversion_md import *
from hydraulic_md import *
from CoolProp.CoolProp import *
from design_md import *
from heat_pump_md import *
from finance_md import *

def Calcul_ccv(kc,kp,Tfo_ch,Tfo_cl,conf,rb,tube_loop):
    g = 9.81
    Patm = 101.325*1000.0
    VFD = True
    if case == -1:
        fluide = 'Water'
        fact_dp = 1.
        fact_cap_ch = 1.
        fact_cap_cl = 1.
        pour = 0
    elif case in [0,1,2,3,5,6,7]:
        fluide  = 'INCOMP::APG-20%'   # propylene - glycol 20
        fact_dp = 1.27
        fact_cap_ch = 0.913
        fact_cap_cl = 0.969
        pour = 20
    elif case == 4:
        fluide  = 'INCOMP::APG-30%'   # propylene - glycol 20
        fact_dp = 1.433
        fact_cap_ch = 0.854
        fact_cap_cl = 0.95
        pour = 30
    Trefk = 277
    muf = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
    rhof = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
    Cpf = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
    Prf = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
    kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
    rhow = PropsSI('D', 'T', Trefk, 'P', Patm, 'Water')
    muw = PropsSI('viscosity', 'T', Trefk, 'P', Patm, 'Water')
    Sg = rhof/rhow
    nuf = muf/rhof     # viscosité cinématique
    nuw = muw/rhow     # viscosité cinématique

    #
    # donnees du problemes
    #
    #
    # parametre du sol
    #
    #
    # parametre du sol
    #
    alj = 0.1 # m2/jr
    alhr = alj/24.0
    ks = 2.75
    To = 10.0
    alsol = alhr/3600 # en m2/s
    my_ground = ground(ksoil=ks,alsoil=alsol,To = To)
    #

    #
    # parametres associés à la configuration des puits
    #
    # diametre des tuyaux dans les boucles pour les 5 champs individuels
    #
    nzones = 3

    tube_head = 3
    #
    # building
    #
    epsilon = 0

    # heat pump
    tube_hoze = [1.25,1.25]
    tube_pac = 1.5
    dist_building = 30
    L_buildingv = np.array([10,66,6,10,6,6,66,6,10,6,2*dist_building])
    tube_buildingv = np.array([3,2,2,1.5,1.25,1.5,2,2,2,3,tube_head])
    gpmv = np.array([45,32,24,16,8,13,26,34,42,50,58])
    #
    #
    # puits
    #

    xc= rb/3.0

    SDR = 11
    di,do = sdr_pipe(tube_loop,SDR)     # choix des tuyaux SDR-11 1.0 po nominal
    extra = 0
    if conf == 'A':
        xc = do/2
    elif conf == 'B':
        xc = rb/3
    elif conf == 'C':
        extra = 5000
        xc = rb - do/2
    #
    # champ de puits
    #
    #
    # champ de puits
    #
    nx = 3 # On suppose que tous les champs ont un seul puits
    ny = 4 # On suppose que tous les champs ont un seul puits
    nb = nx*ny  # nombre de puits pour chaque champ
    d = 6.0     # espacement entre les puits si plus que 1
    nbi = int(nb/ny)     # nombre de trajets parallèles
    # interconnexions
    tube_field = [2,1.5,1.25,tube_loop]
    L_field = [d,d,d,d]
    tube_retour = 2.0
    L_retour = 4*d

    #
    # economic data
    #
    down_payment = 0.20
    t_discount = 0.06                # taux d'actualisation
    inte = 0.04             # interet
    Nccv =  20                  # nombre d,annéees  pour le calcul de la VAN
    Npret = 10              #  number of mortgage loan
    inf_comb = 0.04         # inflation sur combustible
    if kc == 1:
        grout_cost = 85.0            # $ /m3
    elif kc == 1.7:
        grout_cost = 150
    Glycol_cost_litre = 12 # $/l
    kWh_cost  = 0.08
    Drilling_cost = 40  # $/m
    Excavation_cost  = 65 # $/m3
    cout_pompe = 2000
    Et = 1.137
    if kp == 0.4:
        Ct = 248.53
    elif kp == 0.75:
        Ct = 2*248.53
    Cw = 272.70
    Ew = 0.84
    Ci = 50
    Ei = 0.52
    #
    # parametres de conception
    #
    n_annees = 20
    nbloc = 4.0
    ta = n_annees*8760.0
    tm = ta + 730.0
    tf = tm + nbloc
    #
    #
    # loads
    #
    #
    fichier = '..\\data\\loads.txt'
    donnees = np.loadtxt(fichier)
    Chauffage_zone1 = donnees[:,0]
    Climatisation_zone1 = donnees[:,1]
    Chauffage_zone2 = donnees[:,2]
    Climatisation_zone2 = donnees[:,3]
    Chauffage_zone3 = donnees[:,4]
    Climatisation_zone3 = donnees[:,5]
    charge_ch = np.zeros((8760,3))
    charge_ch[:,0] = Chauffage_zone1
    charge_ch[:,1] = Chauffage_zone2
    charge_ch[:,2] = Chauffage_zone3
    charge_cl = np.zeros((8760,3))
    charge_cl[:,0] = Climatisation_zone1
    charge_cl[:,1] = Climatisation_zone2
    charge_cl[:,2] = Climatisation_zone3
    MAX_ch_zone1 = max(Chauffage_zone1)
    MAX_ch_zone2 = max(Chauffage_zone2)
    MAX_ch_zone3 = max(Chauffage_zone3)
    MAX_MBTU_ch_zone1 = BTUhr_W(MAX_ch_zone1)
    MAX_MBTU_ch_zone2 = BTUhr_W(MAX_ch_zone2)
    MAX_MBTU_ch_zone3 = BTUhr_W(MAX_ch_zone3)
    MAX_cl_zone1 = max(Climatisation_zone1)
    MAX_cl_zone2 = max(Climatisation_zone2)
    MAX_cl_zone3 = max(Climatisation_zone3)
    MAX_MBTU_cl_zone1 = BTUhr_W(MAX_cl_zone1)
    MAX_MBTU_cl_zone2 = BTUhr_W(MAX_cl_zone2)
    MAX_MBTU_cl_zone3 = BTUhr_W(MAX_cl_zone3)
    #
    # choice of heat pumps
    #
    # zone 1 : 2 NV060
    #
    EWT_SH = F_C(Tfo_ch) # arbitrary choice minimum temperature
    EWT_SC = F_C(Tfo_cl) # arbitrary choice minimum temperature
    name_zone1 = '..\\data\\wa_nv060.xls'
    df1 = pd.read_excel(name_zone1)
    gpm_zone1,gpm_min_zone1,gpm_max_zone1 = HP_flow_rates(df1)
    cfm_loadv  = df1['CFM_H'].values
    cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
    cfm_load1 = max(cfm_loadv)
    CAP_ch_zone1,Power1,COP_ch_zone1 = WA_heat_pump(EWT_SH,gpm_zone1,cfm_load1,'heating',df1)
    Dppsi1 = fact_dp*WA_hp_psi(EWT_SH,gpm_zone1,df1)
    cfm_loadv  = df1['CFM_C'].values
    cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
    cfm_load1 = max(cfm_loadv)
    CAP_cl_zone1,SC,Power1,EER_zone1 = WA_heat_pump(EWT_SC,gpm_zone1,cfm_load1,'cooling',df1)
    COP_cl_zone1 = W_BTUhr(EER_zone1)
    DpkPa1 = kPa_psi(Dppsi1)
    npac1 = 2


    HP1 = heat_pump(flowrate=m3s_gpm(gpm_zone1),Dp = 1000*DpkPa1,\
        heating_capacity = 1000*W_BTUhr(CAP_ch_zone1), \
        cooling_total_capacity = 1000*W_BTUhr(CAP_cl_zone1),\
        cooling_cop = COP_cl_zone1,heating_cop =  COP_ch_zone1)

    #
    # zone 2 : 2 NV036
    #
    name_zone2 = '..\\data\\wa_nv036.xls'
    df2 = pd.read_excel(name_zone2)
    gpm_zone2,gpm_min_zone2,gpm_max_zone2 = HP_flow_rates(df2)
    cfm_loadv  = df2['CFM_H'].values
    cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
    cfm_load2 = max(cfm_loadv)
    CAP_ch_zone2,Power2,COP_ch_zone2 = WA_heat_pump(EWT_SH,gpm_zone2,cfm_load2,'heating',df2)
    Dppsi2 = fact_dp*WA_hp_psi(EWT_SH,gpm_zone2,df2)
    cfm_loadv  = df2['CFM_C'].values
    cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
    cfm_load2 = max(cfm_loadv)
    CAP_cl_zone2,SC,Power2,EER_zone2 = WA_heat_pump(EWT_SC,gpm_zone2,cfm_load2,'cooling',df2)
    COP_cl_zone2 = W_BTUhr(EER_zone2)
    DpkPa2 = kPa_psi(Dppsi2)
    npac2 = 2

    HP2 = heat_pump(flowrate=m3s_gpm(gpm_zone2),Dp = 1000*DpkPa2,\
        heating_capacity = 1000*W_BTUhr(CAP_ch_zone2), \
        cooling_total_capacity = 1000*W_BTUhr(CAP_cl_zone2),\
        cooling_cop = COP_cl_zone2,heating_cop =  COP_ch_zone2)
    #
    # zone 3 : 2 NV036
    #
    #
    gpm_zone3 = gpm_zone2
    npac3 = npac2
    CAP_cl_zone3 = CAP_cl_zone2
    CAP_ch_zone3 = CAP_ch_zone2
    COP_cl_zone3 = COP_cl_zone2
    COP_ch_zone3 = COP_ch_zone2
    Dppsi3 = Dppsi2
    DpkPa3 = DpkPa2
    #
    # list os the 6 heat pumps for the 3 zones
    #
    HPv = [HP1,HP1,HP2,HP2,HP2,HP2]
    #
    # calcul des capacités e
    #
    CAP_cl = fact_cap_cl*np.array([npac1*CAP_cl_zone1,npac2*CAP_cl_zone2,npac3*CAP_cl_zone3])
    CAP_ch = fact_cap_ch*np.array([npac1*CAP_ch_zone1,npac2*CAP_ch_zone2,npac3*CAP_ch_zone3])
    COPch = np.array([COP_ch_zone1,COP_ch_zone2,COP_ch_zone3])
    COPcl = np.array([COP_cl_zone1,COP_cl_zone2,COP_cl_zone3])
    CAP_cl_kW = W_BTUhr(CAP_cl)
    CAP_ch_kW = W_BTUhr(CAP_ch)
    CAP = sum(np.maximum(CAP_ch_kW,CAP_cl_kW))
    debit_zone = np.array([npac1*HP1.flowrate,npac2*HP2.flowrate,npac3*HP2.flowrate]) # kW
    debit_total = np.sum(debit_zone)

    #

    #
    #
    #  temps d'operation de la pompe de circulation
    #
    vmin = 0.3
    gammav = np.zeros(8760)
    nhv = np.zeros((8760,3))
    for i in range(0,8760):
        gam = 0
        for j in range(0,nzones):
            if charge_ch[i,j] > 0:
    #            gam = charge_ch[i,j]/CAP_ch_kW[j]
    #            if gam > 1:
    #                nhv[i,j] = 1
    #                gam = 1
                gam  = gam + min(charge_ch[i,j]/CAP_ch_kW[j],1)*debit_zone[j]
    #            gam = gam*debit_zone[j]
            else:
    #            gam = charge_cl[i,j]/CAP_cl_kW[j]
    #            if gam > 1:
    #                nhv[i,j] = 1
    #                gam = 1
    #            gam = gam*debit_zone[j]
                gam  = gam + min(charge_cl[i,j]/CAP_cl_kW[j],1)*debit_zone[j]
        gammav[i] = max(gam/debit_total,vmin)
    gammav = np.zeros(8760)
    for i in range(0,8760):
        gam = 0
        for j in range(0,nzones):
            if charge_ch[i,j] > 0:
                x = charge_ch[i,j]/CAP_ch_kW[j]
                if x > 1:
                    nhv[i,j] = 1
                    x = 1
                gam = gam + x*debit_zone[j]
            else:
                x = charge_cl[i,j]/CAP_cl_kW[j]
                if x > 1:
                    nhv[i,j] = 1
                    x = 1
                gam = gam + x*debit_zone[j]
        gammav[i] = max(gam/debit_total,vmin)
    #
    #
    #
    #

    # Calcul des charges au sol horaires en Watts en suposant les COP fixes
    #
    charge_sol_ch = charge_ch*(COPch-1)/COPch*1000.0
    #
    # On assigne le signe moins pour les charges de climatisation
    #
    charge_sol_cl = -charge_cl*(COPcl+1)/COPcl*1000.0
    #
    # q_sol représente les charges horaires côté sol pour les 5 azones mis ensemble
    #
    q_sol = charge_sol_ch + charge_sol_cl
    #
    # calcul des pulses selon la convention de l'ASHRAE pour les 5 systèmes
    #
    charge_totales_ch =  np.sum(charge_sol_ch,axis = 1)
    charge_totales_cl =  np.sum(charge_sol_cl,axis = 1)
    #
    # q_sol représente les charges horaires côté sol pour les 5 azones mis ensemble
    #
    q_sol = charge_totales_ch + charge_totales_cl
    qa,qm_ch,qm_cl,qh_ch,qh_cl =  Pulses_ashrae(q_sol,nbloc = nbloc)

    #
    # ON calcule la puissance au compressur en kWatts
    #
    Wcomp_ch = charge_ch/COPch
    Wcomp_cl = charge_cl/COPcl
    #
    # fin du calcul des charges pour l'application de la methode de l'ASHRAE
    #
    # debut du dimensionnement
    #
    # calcul des débits
    #
    mp = debit_total*rhof
    #
    # calcul des résistances de puits par la ligne source
    #
    ro = do/2.0
    ri = di/2.0
    Rp,Rcond,Rconv = Rp_fct(mp/nb,ro,ri,kp,Trefk,fluide)
    Rpc,Rpac = Rb_linesource(kc,ks,rb,ro,xc)
    Rpb = Rpc + Rp/2.0
    Rpa = Rpac + 2*Rp
    #
    # Calcul des longueurs pour chacun des champs
    #
    CCf = mp*Cpf
    my_bore = borehole(nx=nx,ny=ny,rb = rb,dist = d,Rb=Rpb,CCf=CCf,Rint = Rpa)
    param_conception_ashrae = ashrae_params(qa=qa,qm_heating = qm_ch, qm_cooling=qm_cl,\
            qh_heating = qh_ch,qh_cooling = qh_cl,Tfo_heating =Tfo_ch,\
            Tfo_cooling  = Tfo_cl,n_years = n_annees, n_bloc = nbloc,flag_Tp = 'ASHRAE',flag_inter = False)
    my_field = borefield(params = param_conception_ashrae,ground = my_ground,borehole = my_bore)
    Longueur =  my_field.Compute_L_ashrae()
    Tpv =  my_field.Tp
    Rbsv =  my_field.Rbs
    Lch =  my_field.L_he
    Lcl =  my_field.L_co
    Tfi_ch =  my_field.Tfi_he

#    print('Ttotal length ',Longueur, ' m')
    H = Longueur/nb
    #
    # heaat pumps head calculations

    #
    # hoze kits
    #
    Cv1 = Cv_valve('Hoze kit',tube_hoze[0])
    Cv2 = Cv_valve('Hoze kit',tube_hoze[1])
    qmhr1 = HP1.flowrate*3600  # m3s/hr
    qmhr2 = HP2.flowrate*3600  # m3s/hr
    DPPa1 = Sg*1e5*(qmhr1/Cv1)**2
    DPPa2 = Sg*1e5*(qmhr2/Cv2)**2
    h_hoze1 = DPPa1/(rhof*g) # m de fluide
    h_hoze2 = DPPa2/(rhof*g) # m de fluide
    #
    # strainers
    #
    Cv_str1 = Cv_valve('Y-strainer',tube_hoze[0])
    Cv_str2 = Cv_valve('Y-strainer',tube_hoze[1])
    DPPas1 = Sg*1e5*(qmhr1/Cv_str1)**2
    DPPas2 = Sg*1e5*(qmhr2/Cv_str2)**2
    h_str1 = DPPas1/(rhof*g) # m de fluide
    h_str2 = DPPas2/(rhof*g) # m de fluide
    hhoze = np.array([h_hoze1,h_hoze1,h_hoze2,h_hoze2,h_hoze2,h_hoze2])
    hstr = np.array([h_str1,h_str1,h_str2,h_str2,h_str2,h_str2])
    #
    #  ball valves
    #
    Cv_val1 = Cv_valve('Ball valves',tube_hoze[0])
    Cv_val2 = Cv_valve('Ball valves',tube_hoze[1])
    DPPas1 = Sg*1e5*(qmhr1/Cv_val1)**2
    DPPas2 = Sg*1e5*(qmhr2/Cv_val2)**2
    h_val1 = DPPas1/(rhof*g) # m de fluide
    h_val2 = DPPas2/(rhof*g) # m de fluide
    hvalve = np.array([h_val1,h_val1,h_val2,h_val2,h_val2,h_val2])
    h_hpv = np.zeros(6)
    di,do = sdr_pipe(tube_pac,SDR)
    A1 = pi*di**2/4.0
    Lpipe = 4
    Leqv_hp = np.zeros(6)
    for i in range(0,6):
        q = HPv[i].flowrate
        u = q/A1
        Re = di*u/nuf
        ed = epsilon/di
        f = Colebrook(Re,ed)
        h_pipe = f*((Lpipe)/di)*u**2/(2*g)
        Dp = HPv[i].Dp   # Pascal
        h_hp = Dp/(rhof*g)
        Leqv_hp[i] = h_pipe + fact_dp*(hhoze[i] + hstr[i] + 2*hvalve[i])
        h_hpv[i] = h_hp + Leqv_hp[i]


    Le_string = [['Butt tee-straight'],['Butt 90','Butt 90','Butt reducer'],\
                ['Butt tee-straight'],['Butt 90','Butt reducer'],['Butt reducer','Butt 90'],\
                ['Butt 90'],['Butt tee-straight','Butt 90'],\
                ['Butt tee-straight'],['Butt tee-straight','Butt 90'],\
                ['Butt tee-straight'],['Butt tee-straight','Butt 90']]
    Leqv = np.zeros(11)
    h_pipev = np.zeros(11)
    for i in range(0,11):
        di,do = sdr_pipe(tube_buildingv[i],SDR)
        A1 = pi*di**2/4.0
        debiti = m3s_gpm(gpmv[i])
        u1 = debiti/A1
        Re = di*u1/nuf
        ed = epsilon/di
        f1 = Colebrook(Re,ed)
        Le = 0
        for Le_s in Le_string[i]:
            Leq = sing_head_loss(Le_s,tube_buildingv[i])
            Le = Le + Leq
        Leqv[i] = Le
        h_pipev[i] = f1*((L_buildingv[i]+Le)/di)*u1**2/(2*g) # Dp en metres

    hb1 = sum(h_pipev[0:5]) + h_hpv[5]
    hb2 = sum(h_pipev[0:4])  + h_pipev[9] + h_hpv[4]
    hb3 = sum(h_pipev[0:3])  + sum(h_pipev[8:10]) + h_hpv[3]
    hb4 = sum(h_pipev[0:2])  + sum(h_pipev[7:10]) + h_hpv[2]
    hb5 = h_pipev[0]  + sum(h_pipev[6:10]) + h_hpv[1]
    hb6 =  sum(h_pipev[5:10]) + h_hpv[0]
    #print ('Path first heat pump= ',hb6)
    #print ('Path first heat pump= ',hb5)
    #print ('Path first heat pump= ',hb4)
    #print ('Path second heat pump=',hb3)
    #print ('Path third heat pump= ',hb2)
    #print ('Path last heat pump= ',hb1)
    h_building = max([hb1,hb2,hb3,hb4,hb5,hb6])
    h_header = h_pipev[10]
    #
    #

    #
    # 4) pertes de charges dans la boucle
    #
    debit_loop = debit_total/nb
    Leq_coude = sing_head_loss('Butt 90',tube_loop)
    Leq_ubend = sing_head_loss('Unicoil',tube_loop)
    Leq_convergent = sing_head_loss('Butt tee-branch',tube_loop)      # en pied
    D_loop =  sdr_pipe(tube_loop,SDR)[0]
    A_loop = pi*D_loop**2/4
    u_loop = debit_loop/A_loop
    Re_loop = D_loop*u_loop/nuf
    print('Reynolds loop = ',Re_loop)
    f_loop = Colebrook(Re_loop,0.0)
    L_loop = 2*Longueur/nb
    Lt_loop = L_loop + Leq_ubend +  Leq_convergent + Leq_coude
    h_loop = f_loop*(Lt_loop/D_loop)*u_loop**2/(2*g)  # % Dp en metres
    nbi = nx
    htotal_loop = h_loop
    debiti = debit_total/nx
    for i in range(0,nbi+1):
        i1 = i+1.0
        x = i1/nb
        D =  sdr_pipe(tube_field[i],SDR)[0]
        A = pi*D**2/4.0
        u = debiti/A
        Re = D*u/nuf
        f = Colebrook(Re,0.0)
        h = f*((L_field [i])/D)*u**2/(2*g) #  en metres
        debiti = debiti - debit_loop
        htotal_loop = htotal_loop + h
        # Somme des pertes de charges
        #
    debit_retour = debit_total/nx
    D_retour =  sdr_pipe(tube_retour,SDR)[0]
    A_retour = pi*D_retour**2/4
    u_retour = debit_retour/A_retour
    Re_retour = D_retour*u_retour/nuf
    f_retour = Colebrook(Re_retour,0.0)
    h_retour = f_retour*(L_retour/D_retour)*u_retour**2/(2*g) #  en metres
    htotal_loop = htotal_loop + h_retour
    h_tot = h_building + h_header + htotal_loop
    #
    # choix des pompes
    #  zone 1 E-90 rend w-w 53 %
    #
    rend = 0.53
    W_pompe= mp*g*h_tot
    P_elec = W_pompe/rend
    #
    # facteur de qualité Watts de pompage / kWatts installé
    #
    ratio = P_elec/CAP
    #print (' W /kW ratio is = ',ratio)
    #
    # Consommation de pompage pour chaque zone en W-h
    # Consommation du compresseur  pour chaque zone en W-h
    #
    Consommation_pompe_continu = P_elec*8760.0  # W-hr
    Wchauff = np.sum(Wcomp_ch,axis=0)*1000.0  # Les données de départ étaient en kWatts
    Wclim = np.sum(Wcomp_cl,axis=0)*1000.0
    Wcomp = Wchauff+Wclim
    #
    Energie_pompe1 = np.sum(Consommation_pompe_continu)
    #print ('pumping energy if running all the time  = ', Energie_pompe1/1000, ' kW-hr')
    Energie_totale1 = Energie_pompe1+sum(Wcomp)
    #print ('Total energy if pumps running all the time   =  ' , Energie_totale1/1000, ' kW-hr')
    ratio1 = Energie_pompe1/Energie_totale1
    #print ('Ratio of pumping energy if pumps run all the time  =  ', ratio1)
    #
    # evaluation of the VFD efficiency
    #
    HP = 1.5
    eta_vfd = VFD_efficiency(HP,1)
    Welec = P_elec/eta_vfd
    X = h_building/h_tot
    Z =  fct_lemire(1,X)
    Wshaft_nom = Welec/Z
    Zv = fct_lemire(gammav,X)
    Consommation_pompe_intermitent = sum(Zv*Wshaft_nom)
    Cons2 = 0
    eta_vfdv = np.zeros(8760)
    for i in range(0,8760):
        hi = h_tot*((1-X)*gammav[i]**2 + X)
        eta_vfdv[i] = VFD_efficiency(HP,gammav[i])
        Wi= gammav[i]*mp*g*hi
        Pii = Wi/rend/eta_vfdv[i]
        Cons2 = Cons2 + Pii
    Energie_pompe2 = np.sum(Consommation_pompe_intermitent)
    #print ('pumping energy with VFD   (using Lemire function) =', Energie_pompe2/1000, ' kW-hr')
    Energie_totale2 = Energie_pompe2 + np.sum(Wcomp)
    #print ('Total energy with VFD   (using Lemire function) = ' , Energie_totale2/1000, ' kW-hr')
    ratio2 = Energie_pompe2/Energie_totale2
    #print ('Ratio of pumping energy with VFD (using Lemire function)  = ', ratio2)
    #print ('Pumping energy with VFD estimating VFD efficiency  = ',Cons2/1000, ' kW-hr')
    Energie_totale3 = Cons2 + np.sum(Wcomp)
    #print ('Total energy with VFD estimating VFD efficiency  = ' , Energie_totale3/1000, ' kW-hr')
    ratio3 = Cons2/Energie_totale3
    #print ('ratio of pumping energy with VFD estimating VFD efficiency', ratio3)
    # Calcul des consommation et de l'appel de puissance mensuelles
    #
    # economique analysis
    #
    if VFD:
        Consommation = Energie_totale3/1000 # en kWh
        cout_energie =  kWh_cost*Consommation
    else:
        Consommation = Energie_totale1/1000 # en kWh
        cout_energie =  kWh_cost*Consommation
    #
    #

    C_pac = 1949.5*CAP**0.665
    # cout tuyau + antigel
    #

    # boucle
    di = sdr_pipe(tube_loop,SDR)[0]
    A_loop = pi*di**2/4
    C_loop = Ct*m_ft(tube_loop/12)**Et
    L_loop = 2*Longueur
    V_loop = A_loop*L_loop
    C_tuyau =  L_loop*C_loop
    C_tuyaub = C_tuyau
    Volume_antigel =  V_loop
    #
    # header
    C_conn = Ct*m_ft(tube_head/12)**Et
    dj = sdr_pipe(tube_head,SDR)[0]
    Aj = pi*dj**2/4
    L_head = L_buildingv[len(L_buildingv)-1]
    Vj = Aj*L_head
    cout_conn = L_head*C_conn
    V_conn = Vj
    C_iso = Ci*m_ft(tube_head/12)**Ei
    cout_iso = L_head*C_iso
    cout_welding = 0
    # field connexions
    nc = len(L_field)
    for i in range(0,nc):
        C_conn = Ct*m_ft(tube_field[i]/12)**Et
        dj = sdr_pipe(tube_field[i],SDR)[0]
        Aj = pi*dj**2/4
        Vj = Aj*L_field[i]
        C_iso = Ci*m_ft(tube_field[i]/12)**Ei
        cout_iso = cout_iso + nbi*L_field[i]*C_iso*2
        cout_conn = cout_conn + nbi*L_field[i]*C_conn*2
        V_conn = V_conn + nbi*Vj*2
    for i in range(0,nc-1):
        C_weld = Cw*m_ft(tube_field[i]/12)**Ew
        cout_welding = cout_welding + nbi*2*C_weld
    C_weld = Cw*m_ft(tube_head/12)**Ew
    cout_welding = cout_welding + 2*(nbi-1)*C_weld
    D_retour =  sdr_pipe(tube_retour,SDR)[0]
    C_conn = Ct*m_ft(tube_retour/12)**Et
    C_iso = Ci*m_ft(tube_retour/12)**Ei
    A_retour = pi*D_retour**2/4
    V_retour = A_retour*L_retour
    V_conn = V_conn + nbi*V_retour
    cout_iso = cout_iso + nbi*L_retour*C_iso
    cout_conn = cout_conn + nbi*L_retour*C_conn
    C_tuyau = C_tuyau + cout_conn + cout_welding + cout_iso
    C_tuyaub = C_tuyaub + cout_conn
    Cost_weld = cout_welding
    Cost_iso = cout_iso
    Volume_antigel = Volume_antigel + V_conn
    #
    # building connextion
    nc = len(L_buildingv)
    cout_build = 0
    cout_iso_build = 0
    V_build  = 0
    cout_welding_build = 0
    for i in range(0,nc):
        C_build = Ct*m_ft(tube_buildingv[i]/12)**Et
        C_iso = Ci*m_ft(tube_buildingv[i]/12)**Ei
        dj = sdr_pipe(tube_buildingv[i],SDR)[0]
        Aj = pi*dj**2/4
        Vj = Aj*L_buildingv[i]
        cout_build = cout_build + L_buildingv[i]*C_build
        cout_iso_build = cout_iso_build + L_buildingv[i]*C_iso
        V_build = V_build + Vj
    for i in range(0,nc-1):
        C_weld = Cw*m_ft(tube_buildingv[i]/12)**Ew
        cout_welding_build = cout_welding_build + C_weld
    Cost_weld = Cost_weld +  cout_welding_build
    Cost_iso = Cost_iso + cout_iso_build
    C_tuyau = C_tuyau + cout_build +   cout_iso_build + cout_welding_build
    C_tuyaub = C_tuyaub + cout_build
    Volume_antigel = Volume_antigel + V_build
    #
    # cout de l'antigel
    #
    C_antigel = Glycol_cost_litre*Volume_antigel*1000*pour/100
    #
    # cout excavation
    #

    Lexc = nbi*sum(L_field) + dist_building
    A_trench = 1*2
    Vexc = Lexc*A_trench
    C_forage = Drilling_cost*Longueur
    C_exc =  Excavation_cost*Vexc
    Ac = pi*rb**2 - 2*pi*ro**2
    Vc = Longueur*Ac
    C_coulis = grout_cost*Vc
    C_pompe = cout_pompe
    if VFD:
        C_vfd = 1000
    else:
        C_vfd = 0
    Cinit = extra + C_pac  + C_forage   + C_exc + C_tuyau + C_coulis + C_antigel + C_pompe + C_vfd
    #print('cout initial = ',Cinit)
    #print('pompe = ',C_pompe)
    #print('pac = ',C_pac)
    #print('tuyau = ',C_tuyau)
    #print('antigel = ',C_antigel)
    #print('excr = ',C_exc)
    #print('coulis = ',C_coulis)
    #print('forage = ',C_forage)
    Invest = Cinit       # cout total
    dpa = down_payment*Invest            # paiement initial
    mi = Invest - dpa
    paiementi = mi/pwf(Npret,0,inte)     # paiment annuel
    paiement = paiementi
    #
    CCV = dpa + paiementi*pwf(Npret,0,t_discount)+cout_energie*pwf(Nccv,inf_comb,t_discount)
    print ('CCV = ',CCV)
    return Cinit,CCV

input_file = '..\\data\\Exemple12_18_input.xlsx'
df_input = pd.read_excel(input_file)
for case in range(0,8):
    kc = df_input['kgrout'][case]
    kp = df_input['kplastic'][case]
    Tfo_ch = df_input['EWT_min'][case]
    Tfo_cl = df_input['EWT_max'][case]
    conf  = df_input['conf'][case]
    rb = df_input['radius'][case]
    tube_loop = df_input['loop_diameter'][case]
    Cinit,ccv = Calcul_ccv(kc,kp,Tfo_ch,Tfo_cl,conf,rb,tube_loop)
    print('case ',case)
    print('Initial cost :',Cinit,'CCV :',ccv)