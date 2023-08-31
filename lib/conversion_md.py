import numpy as np
#
#
#
#
def help_conversion():
    print('flow')
    print ('ls_gpm (gpm=1):')
    print ('mcs_gpm (gpm=1):')
    print ('gpm_mcs (mcs=1):')
    print ('gpm_ls (ls=1):')
    print('temperature)')
    print ('F_C (C=1):')
    print ('C_F (F=1):')
    print ('ls_kW_gpm_tonne (gpm_tonne): # ls_kW_gpm_tonne : transforme des gpm par tonne en l/s par kWattts')
    print ('kWh_MTOE (mtoe = 1):')
    print ('MJ_MTOE (mtoe = 1):')
    print ('W_m2K_BTU_hrft2F (BTU_hrft2F=1):')
    print ('W_mK_BTU_hrftF (BTU_hrftF=1):')
    print ('W_BTUhr (BTU_hr=1):')
    print ('W_kW_hp_tonne (hp_tonne):')
    print ('J_BTU (btu =1):')
    print ('MJ_m2_kBTU_f2 (kbtu_f2 =1):')
    print ('W_hp (hp = 1):')
    print('\t pressure')
    print ('Pa_inw (inw = 1):')
    print ('inw_Pa (Pa = 1):')
    print ('Pa_ftw (ftw = 1):')
    print ('mw_psi (psi = 1):')
    print ('kPa_psi (psi = 1):')
    print ('Pa_bar(bar = 1):')
    print('psi_bar(bar = 1):')
    print('def psi_kPa(kpa = 1')
    print('bar_psi(psi = 1):')
    print('\t distance-area-volume')
    print ('m_ft (ft =1):')
    print ('m2_ft2(ft =1):')
    print ('ft_m (m =1):')
    print ('ft2_m2 (m =1):')
    print ('kg_lb(lb = 1):')
    print ('lb_kg(kg = 1):')
#
    print('\t  conductivity - resistance')
    print ('W_m2K_BTU_hrft2F (U=1):')      # U - factor
    print ('BTU_hrft2F_W_m2K(U=1):')      # U - factor
    print ('Rsi_R(R =1):')      # R - factor
    print ('R_Rsi(Rsi=1):')      # R - factor
    print ('W_mK_BTU_hrftF(kip=1):')           # conductivity
    print ('BTU_hrftF_W_mK(ksi = 1):')
    print('\t  # mass heat capacity')
    print ('BTU_lbF_J_kgK(Cp = 1):')
    print ('J_kgK_BTU_lbF(Cp = 1):')




#######
# flow
def m3s_cfm (cfm=1):
    ls = cfm*ls_cfm()
    mcs = ls/1000.0
    return mcs
def m3s_gpm (gpm=1):
    ls = gpm*ls_gpm()
    mcs = ls/1000.0
    return mcs
def cfm_m3s (mcs=1):
    ls = 1000.0*mcs
    cfm = cfm_ls(ls)
    return cfm
def gpm_m3s (mcs=1):
    ls = 1000.0*mcs
    gpm = gpm_ls(ls)
    return gpm
def gpm_ls (ls=1):
    gpm = ls/0.0630902
    return gpm
def ls_gpm (gpm=1):
    ls = gpm*0.0630902
    return ls
def ls_cfm (cfm=1):
    ls = cfm*0.4719474
    return ls
def cfm_ls (ls=1):
    cfm = ls/0.4719474
    return cfm
# temperature
#
def F_C (C=1):
    F = C*1.8 + 32.0
    return F
def C_F (F=1):
    C = (F-32.0)/1.8
    return C
#
#
# varia
#
def ls_kW_gpm_ton(gpm_ton): # ls_kW_gpm_tonne : transforme des gpm/ton in  l/s per kW
    ls_ton = gpm_ton*ls_gpm()
    kW_ton = W_BTUhr()*12.0
    ls_kW = ls_ton/kW_ton
    return ls_kW

# POwer and Energy
def W_HP (HP=1):
    W = HP*745.6999
    return W
def HP_W (W=1):
    HP = W/745.6999
    return HP
def W_BTUhr (BTU_hr=1):
    W = BTU_hr*0.2928104
    return W
def BTUhr_W (W = 1):
    BTU_hr = W/0.2928104
    return BTU_hr
def W_kW_hp_ton(hp_ton):
    W_ton = hp_ton*W_hp()
    W_MBTUhr = W_ton/12.0
    W_kW = W_MBTUhr/W_BTUhr()
    return W_kW
def J_Wh(Wh =1):
    joule = Wh*3600.0
    return joule
def Wh_J(J =1):
    Wh = J/3600.0
    return Wh
def J_BTU (btu =1):
    joule = btu*1055.056
    return joule
def BTU_J (J =1):
    btu = J/1055.056
    return btu
def BTU_lb_J_kg(J_kg =1):
    btu = BTU_J(J_kg)/lb_kg()
    return btu
def J_kg_BTU_lb(btu_lb =1):
    btu = J_BTU(btu_lb)*lb_kg()
    return btu
def MJ_m2_kBTU_f2 (kbtu_f2 =1):
    joule_f2 = kbtu_f2*1000*1055.056
    mjoule_m2 = joule_f2/m_ft()/m_ft()/1e6
    return mjoule_m2
def W_hp (hp = 1):
    W = hp*745.6999
    return W
def kWh_MTOE(mtoe = 1):
    kwh = mtoe*1.163*1e10
    return kwh
def MJ_MTOE(mtoe = 1):
    kwh = mtoe*4.1868*1e10
    return kwh
def BBL_MJ(mj = 1):
    bbl = mj*1.64e-4
    return bbl
def MJ_BBL(bbl = 1):
    mj = bbl*6.1e3
    return mj
def BBL_toe(t = 1):
    bbl = t*7.33
    return bbl
def toe_BBL(bbl = 1):
    mj = bbl*0.136
    return mj
def BBL_m3(m3 = 1):
    bbl = m3*6.43e-3
    return bbl
def m3_BBL(bbl = 1):
    m3 = bbl*155
    return m3



# pressure
def Pa_ftw (ftw = 1):
    Pa = ftw *2990.1
    return W
def Pa_inw (inw = 1):
    Pa = inw/12*2990.1
    return Pa
def inw_Pa (Pa = 1):
    inw = Pa/2990.1*12
    return inw
def ftw_psi(psi = 1):
    fw = ft_m(mw_psi(psi))
    return fw
def inw_psi(psi = 1):
    iw = ft_m(mw_psi(psi))*12
    return iw
def mw_psi (psi = 1):
    y = 0.702829*psi
    return y
def kPa_psi (psi = 1):
    y = 6.89457*psi
    return y
def Pa_bar(bar = 1):
    Pa =  bar*1e5
    return Pa
def psi_bar(bar = 1):
    psi =  bar*14.50377
    return psi
def psi_kPa(kpa = 1):
    psi =  kpa*14.50377/100
    return psi
def bar_psi(psi = 1):
    bar =  psi/14.50377
    return bar
def kPa_psi(psi = 1):
    kPa =  psi/.1450377
    return kPa
# distance-area-volume
def m_ft (ft =1):
    m = ft*0.3048
    return m
def m2_ft2(ft =1):
    m = ft*0.3048**2
    return m
def m3_ft3(ft =1):
    m = ft*0.3048**3
    return m
def m3_gal(gal =1):
    m = gal*0.003785412
    return m
def gal_m3(m = 1):
    gal = m/0.003785412
    return gal
def ft_m (m =1):
    ft = m/0.3048
    return ft
def ft2_m2 (m =1):
    ft = m/0.3048**2
    return ft
def ft3_m3 (m =1):
    ft = m/0.3048**3
    return ft
def gal_l (l =1):
    gal = l/3.785412
    return gal
def l_gal (g =1):
    l = g*3.785412
    return l
def kg_lb(lb = 1):
    kg = lb*0.4535924
    return kg
def lb_kg(kg = 1):
    lb = kg/0.4535924
    return lb
#
# vavles Cv
#
def gpm_psi_m3h_bar(kv = 1):
    a = gpm_m3s()/3600
    b = psi_bar()
    cv = kv*a/np.sqrt(b)
    return cv
def gpm_psi_ls_kPa(kv = 1):
    a = gpm_ls()
    b = psi_kPa()
    cv = kv*a/np.sqrt(b)
    return cv
def m3h_bar_gpm_psi(cv = 1):
    a = m3s_gpm()*3600
    b = bar_psi()
    kv = cv*a/np.sqrt(b)
    return kv
def ls_kPa_gpm_psi(cv = 1):
    a = ls_gpm()
    b = kPa_psi()
    kv = cv*a/np.sqrt(b)
    return kv
def m2hr_bar_gpm_psi(cv = 1):
    kv = cv*0.8649894534152124
    return kv
#
# conductivity - resistance




def W_m2K_BTU_hrft2F (U=1):      # U - factor
    W = U*W_BTUhr()*ft_m()*ft_m()*1.8
    return W
def BTU_hrft2F_W_m2K(U=1):      # U - factor
    Wn = U/W_m2K_BTU_hrft2F()
    return Wn

def Rsi_R(R =1):      # R - factor
    Rn = R/W_m2K_BTU_hrft2F()
    return Rn
def R_Rsi(Rsi=1):      # R - factor
    Rn = Rsi*W_m2K_BTU_hrft2F()
    return Rn

def W_mK_BTU_hrftF(kip=1):           # conductivity
    k = kip*W_BTUhr()*ft_m()*1.8
    return k
def BTU_hrftF_W_mK(ksi = 1):
    kn = ksi/W_mK_BTU_hrftF()
    return(kn)

# mass heat capacity
def BTU_lbF_J_kgK(Cp = 1):
    c1 = BTU_J()
    c2 = 2.20462 # lb/kg
    c = c1/c2/1.8
    Cpn = Cp*c
    return(Cpn)

def J_kgK_BTU_lbF(Cp = 1):
    c1 = BTU_J()
    c2 = 2.20462 # lb/kg
    c = c1/c2/1.8
    Cpn = Cp/c
    return(Cpn)


# volumetric heat capacity

def BTU_ft3F_J_m3K(C = 1):
    c1 = BTU_J()
    c2 = m_ft()
    c = c1*c2**3/1.8
    Cn = C*c
    return(Cn)

def J_m3K_BTU_ft3F(C = 1):
    c1 = BTU_J()
    c2 = m_ft()
    c = c1*c2**3/1.8
    Cn = C/c
    return(Cn)

# density
def lb_ft3_kg_m3(rho = 1):
    c1 = 2.20462 # lb/kg
    c2 = m_ft()
    c = c1*c2**3
    rhon = rho*c
    return(rhon)
def kg_m3_lb_ft3(rho = 1):
    c1 = 2.20462 # lb/kg
    c2 = m_ft()
    c = c1*c2**3
    rhon = rho/c
    return(rhon)

# viscosity

def lbf_ft2_s_Pa_s(mu = 1):
    c1 =  0.22480894387096  # lbf/N
    c2 = m_ft()
    c = c1*c2**2
    mun = mu*c
    return(mun)

def Pa_s_lbf_ft2_s(mu = 1):
    c1 =  0.22480894387096  # lbf/N
    c2 = m_ft()
    c = c1*c2**2
    mun = mu/c
    return(mun)

def lbm_ft_hr_Pa_s(mu = 1):
    c =  0.0004133788701  # lbf/N
    mun = mu/c
    return(mun)

def Pa_s_lbm_ft_hr(mu = 1):
    c =  0.0004133788701  # lbf/N
    mun = mu*c
    return(mun)
