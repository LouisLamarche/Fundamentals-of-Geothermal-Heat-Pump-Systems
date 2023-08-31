#coding: utf-8
from geothermal_md import *
import  numpy as np
from conversion_md import *
from hydraulic_md import *
from CoolProp.CoolProp import *
from design_md import *
import pandas as pd
from heat_pump_md import *


jm = np.array([17,47,75,105,135,162,198,228,258,288,318,344])# jour type de chaque mois
nm = np.array([31,28,31,30,31,30,31,31,30,31,30,31])                # nombre de jours par mois
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])  # nombre d'heures dans chaque mois
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]) # nombre d'heures écoulées après chaque mois
g = 9.81
Patm = 101.325*1000.0
cas = 2
if cas == 1:
    fluide = 'Water'
    fact_dp = 1.
    fact_cap_ch = 1.
    fact_cap_cl = 1.
elif cas ==2:
    fluide  = 'INCOMP::APG-20%'   # propylene - glycol 20
    fact_dp = 1.27
    fact_cap_ch = 0.913
    fact_cap_cl = 0.969

Trefk = 277
muf = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
rhof = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cpf = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Prf = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
rho_eau = PropsSI('D', 'T', Trefk, 'P', Patm, 'Water')
mu_eau = PropsSI('viscosity', 'T', Trefk, 'P', Patm, 'Water')
Sg = rhof/rho_eau
nuf = muf/rhof     # viscosité cinématique
nu_eau = mu_eau/rho_eau     # viscosité cinématique
#
# data
#
# soil parameters
#
alj = 0.1 # m2/jr
alhr = alj/24.0
ks = 2.75
To = 10.0
alsol = alhr/3600 # en m2/s
my_ground = ground(ksoil=ks,alsoil=alsol,To = To)
#
#
#
nzones = 3
tube_loopv = np.array([1,1,1])
tube_maiv = np.array([1.25,1.25,1.25])
tube_headv = np.array([2,1.25,1.25])
L_head = 60

#
#
# puits
#
kc = 1.0
rb = 0.15/2.0
xc= rb/3.0
kp = 0.4
SDR = 11
div = np.zeros(nzones)
dov = np.zeros(nzones)
for i in range(0,nzones):
    div[i],dov[i] = sdr_pipe(tube_loopv[i],SDR)     # choix des tuyaux SDR-11 1.0 po nominal
#
# champ de puits
#
nxv = np.array([2,3,3]) # On suppose que tous les champs ont un seul puits
nyv = np.array([2,1,1]) # On suppose que tous les champs ont un seul puits
nbb = nxv*nyv  # nombre de puits pour chaque champ
d = 6.0     # espacement entre les puits si plus que 1
# interconnexions
tube1 = np.array([1.5,1.25,tube_loopv[0]])
tube2 = np.array([1.25,tube_loopv[1]])
tube3 = np.array([1.25,tube_loopv[2]])
tube_field = [tube1,tube2,tube3]
L1 = np.array([d,2*d,d])
L2 = np.array([d,d])
L3 = np.array([d,d])
L_field = [L1,L2,L3]
#
# parametres de conception
#
Tfo_ch = 0.0    # valeur de design
Tfo_cl = 25.0   # valeur de design
n_annees = 20
nbloc = 4.0
ta = n_annees*8760.0
tm = ta + 730.0
tf = tm + nbloc
#
# loads
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
Dppsi1 = WA_hp_psi(EWT_SH,gpm_zone1,df1)
cfm_loadv  = df1['CFM_C'].values
cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
cfm_load1 = max(cfm_loadv)
CAP_cl_zone1,SC,Power1,EER_zone1 = WA_heat_pump(EWT_SC,gpm_zone1,cfm_load1,'cooling',df1)
COP_cl_zone1 = W_BTUhr(EER_zone1)
DpkPa1 = kPa_psi(Dppsi1)
npac1 = 2


HP1 = heat_pump(flowrate=m3s_gpm(gpm_zone1),Dp = 1000*fact_dp*DpkPa1,\
    heating_capacity = 1000*fact_cap_ch*W_BTUhr(CAP_ch_zone1), \
    cooling_total_capacity = 1000*fact_cap_cl*W_BTUhr(CAP_cl_zone1),\
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
Dppsi2 = WA_hp_psi(EWT_SH,gpm_zone2,df2)
cfm_loadv  = df2['CFM_C'].values
cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
cfm_load2 = max(cfm_loadv)
CAP_cl_zone2,SC,Power2,EER_zone2 = WA_heat_pump(EWT_SC,gpm_zone2,cfm_load2,'cooling',df2)
COP_cl_zone2 = W_BTUhr(EER_zone2)
DpkPa2 = kPa_psi(Dppsi2)
npac2 = 2


HP2 = heat_pump(flowrate=m3s_gpm(gpm_zone2),Dp = 1000*fact_dp*DpkPa2,\
    heating_capacity = 1000*fact_cap_ch*W_BTUhr(CAP_ch_zone2), \
    cooling_total_capacity = 1000*fact_cap_cl*W_BTUhr(CAP_cl_zone2),\
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
DpkPa = fact_dp*np.array([DpkPa1,DpkPa2,DpkPa3])
npacv = np.array([npac1,npac2,npac3])
debit_zone = np.array([npac1*HP1.flowrate,npac2*HP2.flowrate,npac3*HP2.flowrate]) # kW
debit_total = np.sum(debit_zone)

#
print('maximum zone1 (Heating) = ',MAX_MBTU_ch_zone1,' MBTU/hr')
print('Capacity zone1 (Heating) = ',CAP_ch[0],' MBTU/hr')
print('maximum zone2 (Heating) = ',MAX_MBTU_ch_zone2,' MBTU/hr')
print('Capacity zone2 (Heating) = ',CAP_ch[1],' MBTU/hr')
print('maximum zone3 (Heating) = ',MAX_MBTU_ch_zone3,' MBTU/hr')
print('Capacity zone3 (Heating) = ',CAP_ch[2],' MBTU/hr')
print('maximum zone1 (Cooling) = ',MAX_MBTU_cl_zone1,' MBTU/hr')
print('Capacity zone1 (Cooling) = ',CAP_cl[0],' MBTU/hr')
print('maimum zone2 (Coo,ing) = ',MAX_MBTU_cl_zone2,' MBTU/hr')
print('Capacity zone2 (Cooling) = ',CAP_cl[1],' MBTU/hr')
print('maximum zone3 (Cooling) = ',MAX_MBTU_cl_zone3,' MBTU/hr')
print('Capacity zone3 (Cooling) = ',CAP_cl[2],' MBTU/hr')

#
#
#  temps d'operation de la pompe de circulation
#  On a fait l'hypothese que la pac choisie a une capacité
#  égale à la capacité maximale nécessaire
#  En pratique on choisie une pac où la capacité est plus grande
#  que la demande max
#
temps_operation = np.zeros((8760,nzones))
for i in range(0,nzones):
    for j in range(0,8760):
        if charge_ch[j,i] > 0:
            temps_operation[j,i] = min(charge_ch[j,i]/CAP_ch_kW[i],1)
        else:
            temps_operation[j,i] = min(charge_cl[j,i]/CAP_cl_kW[i],1)
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
qa = np.zeros(nzones)
qm_ch = np.zeros(nzones)
qm_cl = np.zeros(nzones)
qh_ch = np.zeros(nzones)
qh_cl = np.zeros(nzones)
for i in range(0,nzones):
    qa[i],qm_ch[i],qm_cl[i],qh_ch[i],qh_cl[i] =  Pulses_ashrae(q_sol[:,i],nbloc = nbloc)

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
CAP = np.maximum(CAP_ch_kW,CAP_cl_kW) # maximum capacity en kW ( voir hypothese sur capacite)
mp = debit_zone*rhof
#
# calcul des résistances de puits par la ligne source
#
Rpb = np.zeros(nzones)
Rpa = np.zeros(nzones)
for i in range(0,nzones):
    ro = dov[i]/2.0
    ri = div[i]/2.0
    Rp,Rcond,Rconv = Rp_fct(mp[i]/nbb[i],ro,ri,kp,Trefk,fluide)
    Rpc,Rpac = Rb_linesource(kc,ks,rb,ro,xc)
    Rpb[i] = Rpc + Rp/2.0
    Rpa[i] = Rpac + 2*Rp
#
# Calcul des longueurs pour chacun des champs
#
Longueur = np.zeros(nzones)
Tpv = np.zeros(nzones)
Rbsv = np.zeros(nzones)
for i in range(0,nzones):
    CCf = mp[i]*Cpf
    my_bore = borehole(nx=nxv[i],ny=nyv[i],rb = rb,dist = d,Rb=Rpb[i],CCf=CCf,Rint = Rpa[i])
    param_conception_ashrae = ashrae_params(qa=qa[i],qm_heating = qm_ch[i], qm_cooling=qm_cl[i],\
            qh_heating = qh_ch[i],qh_cooling = qh_cl[i],Tfo_heating =Tfo_ch,\
            Tfo_cooling  = Tfo_cl,n_years = n_annees, n_bloc = nbloc,flag_Tp = 'ASHRAE',flag_inter = False)
    my_field = borefield(params = param_conception_ashrae,ground = my_ground,borehole = my_bore)
    Longueur[i] =  my_field.Compute_L_ashrae()
    Tpv[i] =  my_field.Tp
    Rbsv[i] =  my_field.Rbs
print('la longueur est ',Longueur, ' m')
H = Longueur/nbb
print(' Total length is '+ str(sum(Longueur)) + ' m')
print(' Bore height  are ',H , ' m')


#
# calcul des pertes de charges pour chaque boucle
#
# tableau (diapo 21 cours 6)  pour les Ball valve, on pourrait rajouter les autres si nécesaire
#

#
h_tot = np.zeros(nzones)
h_pacv = DpkPa*1000/(rhof*g)
h_maiv = np.zeros(nzones)
h_headv = np.zeros(nzones)
h_loopv = np.zeros(nzones)
h_hpv = np.zeros(nzones)
for i in range(0,nzones):
    nbore = nbb[i]
    npac = npacv[i]
#
#   Calcul des pertes de charges  PAC
#
    debit_pac = debit_zone[i]/npac
    debit_loop = debit_zone[i]/nbore
#
#   Calcul des pertes de charges entre PAC et collecteur
#


    tube_mai = tube_maiv[i]
    Cv_ball_valve = Cv_valve('Ball valves',tube_mai)   # m3/hr/bar
    Cv_check_valve = Cv_valve('Swing check valves',tube_mai)  # m3/hr/bar
    Cv_strainer =  Cv_valve('Y-strainer',tube_mai) # m3/hr/bar
    Cv_hoze_kit  = Cv_valve('Hoze kit',tube_mai)  # m3/hr/bar
    qmhr = debit_pac*3600  # m3/heure
    DPPa0 = Sg*1e5*(qmhr/Cv_ball_valve)**2
    DPPa1 = Sg*1e5*(qmhr/Cv_check_valve)**2
    DPPa2 = Sg*1e5*(qmhr/Cv_strainer)**2
    DPPa3 = Sg*1e5*(qmhr/Cv_hoze_kit)**2
    h_sing_pac = (2*DPPa0 + DPPa1 + DPPa2+DPPa3)/(rhof*g) # m de fluide
    h_maiv[i] = h_pacv[i] + h_sing_pac*fact_dp

#
# 3) Calcul des pertes de charges dans le header
#
    tube_head = tube_headv[i]
    D_head =  sdr_pipe(tube_head,SDR)[0]
    Leq_coude = sing_head_loss('Butt 90',tube_head)
    A_head = pi*D_head**2/4
    u_head = debit_zone[i]/A_head
    Re_head = D_head*u_head/nuf
    f_head = Colebrook(Re_head,0.0)
    Lt_head = L_head
    h_head = f_head*(Lt_head/D_head)*u_head**2/(2*g)
    h_headv[i] = h_head

#
# 4) pertes de charges dans la boucle
#
    tube_loop = tube_loopv[i]
    Leq_coude = sing_head_loss('Butt 90',tube_loop)
    Leq_ubend = sing_head_loss('Unicoil',tube_loop)
    Leq_convergent = sing_head_loss('Butt tee-branch',tube_loop)      # en pied
    D_loop =  sdr_pipe(tube_loop,SDR)[0]
    A_loop = pi*D_loop**2/4
    u_loop = debit_loop/A_loop
    Re_loop = D_loop*u_loop/nuf
    f_loop = Colebrook(Re_loop,0.0)
    L_loop = 2*Longueur[i]/nbore
    Lt_loop = L_loop + Leq_ubend +  Leq_convergent + Leq_coude
    h_loop = f_loop*(Lt_loop/D_loop)*u_loop**2/(2*g)  # % Dp en metres

    tube = tube_field[i]
    Long = L_field[i]
    htotal_loop = 0
    debit_volumique = debit_zone[i]
    for j in range(0,nbore-1):
        debit_volumique = debit_volumique - debit_loop          #  débit qui passe dans la section i
        D =  sdr_pipe(tube[j],SDR)[0]
        A = pi*D**2/4.0
        u = debit_volumique/A
        Re = D*u/nuf
        f = Colebrook(Re,0.0)
        Leq = 0
        h = f*((Long[j]+Leq)/D)*u**2/(2*g) #  en metres
        htotal_loop = htotal_loop + h
    h_loopv[i]  = h_loop + htotal_loop

    #
    # Somme des pertes de charges
    #
h_tot = h_loopv + h_headv + h_maiv
#
# choix des pompes
#  zone 1 E-90 rend w-w 40.2 %
#  zone 2 E-1535 rend w-w 30 %
#
rend1 = 0.4
rend2 = 0.32
rend3 = 0.32
rend = np.array([rend1,rend2,rend3])
W_pompe= mp*g*h_tot
P_elec = W_pompe/rend
#
# facteur de qualité Watts de pompage / kWatts installé
#
ratio = P_elec/CAP
print ('W /kW ratio are = ')
print (ratio)
#
# Consommation de pompage pour chaque zone en W-h
# Consommation du compresseur  pour chaque zone en W-h
#
Consommation_pompe_continu = P_elec*8760.0  # W-hr
Consommation_pompe_intermittent = np.sum(P_elec*temps_operation,axis=0)  # W-hr
Wchauff = np.sum(Wcomp_ch,axis=0)*1000.0  # Les données de départ étaient en kWatts
Wclim = np.sum(Wcomp_cl,axis=0)*1000.0
Wcomp = Wchauff+Wclim
#
Energie_pump1 = sum(Consommation_pompe_continu)
print ('pumping energy if running all the time  = ', Energie_pump1/1000, ' kW-hr')
Energie_pump2 = sum(Consommation_pompe_intermittent)
print ('pumping energy if stopped with compressor  = ' , Energie_pump2/1000, ' kW-hr')
Energie_totale1 = sum(Consommation_pompe_continu+Wcomp)
print ('Total energy if pumps running all the time   = ' , Energie_totale1/1000, ' kW-hr')
Energie_totale2 = sum(Consommation_pompe_intermittent+Wcomp)
print ('Total energy if pumps are stopped   = ' , Energie_totale2/1000, ' kW-hr')
ratio1 = sum(Consommation_pompe_continu)/Energie_totale1
ratio2 = sum(Consommation_pompe_intermittent)/Energie_totale2
print ('Ratio of pumping energy if pumps run all the time  = ', ratio1)
print ('Ratio of pumping energy if pumps are stopped  = ', ratio2 )
