#coding: utf-8
import numpy as np
from scipy.interpolate import interpn,griddata


def  cherche_index(xi,x):
    """ cherche l'index où x(i) <= xi < x(i+1)"""
    err = 0
    ok = False
    i = 0
    if (xi < x[0] or xi > x[len(x)-1]):
        err =-1
        i = np.nan
    elif xi == x[len(x)-1]:
        i = len(x)-2
    else:
        while not ok:
            if (xi>=x[i]) & (xi<x[i+1]):
                ok = True
            else:
                i=i+1
                if i >= len(x):
                    ok = True
                    err =-1
                    i = nan
    return i



class heat_pump:

    """

    Attributes
    ----------

    """

    def __init__(self, flowrate=0, heating_capacity=0, heating_cop=3, \
        cooling_total_capacity=0, cooling_sensible_capacity=0, \
        cooling_cop = 3,Dp =  0 ):
        self.flowrate = float(flowrate)      # m3/s
        self.heating_capacity = float(heating_capacity)      # Watts
        self.heating_cop = float(heating_cop)
        self.cooling_total_capacity = float(cooling_total_capacity)      # Watts
        self.cooling_sensible_capacity = float(cooling_sensible_capacity)      # Watts
        self.cooling_cop = float(cooling_cop)  #
        self.Dp = float(Dp)      # pressure drop  Pascals



# air_flow_corr
# air_temp_cooling_corr
# air_temp_heating_corr
#
# correction factors for air flow
#
def air_flow_corr(x_air,df,mode = 'cooling'):
    x_air_flow = df['FLOW'].values
    tc_air_flow= df['TC_CORR'].values
    sc_air_flow= df['SC_CORR'].values
    hc_air_flow= df['HC_CORR'].values
    powc_air_flow= df['C_POW_CORR'].values
    powh_air_flow= df['H_POW_CORR'].values
    x_air = x_air*100
    if mode == 'cooling':
        tc_corr = np.interp(x_air,x_air_flow,tc_air_flow)
        sc_corr = np.interp(x_air,x_air_flow,sc_air_flow)
        pow_corr = np.interp(x_air,x_air_flow,powc_air_flow)
        return tc_corr,sc_corr,pow_corr
    elif mode == 'heating':
        hc_corr = np.interp(x_air,x_air_flow,hc_air_flow)
        pow_corr = np.interp(x_air,x_air_flow,powh_air_flow)
        return hc_corr,pow_corr
#
# correction factors for air flow Water_furnace
#


def air_temp_cooling_corr(dry_bulb,wet_bulb,df):
    #
    # correction factors for air temperature cooling
    # dry_bulb and wet_bulb are in F
    #
    x_sens_cool =  df['EA_WBT'].values
    y_sens_cool =  df['EA_DBT'].values
    z_sens_cool =  df['SC_CORR'].values
    tc_corr =  df['TC_CORR'].values
    x_wet_bulb = np.unique(x_sens_cool)
    tc_wet_bulb = np.unique(tc_corr)
    pow_corr =  df['POW_CORR'].values
    power_wet_bulb = np.unique(pow_corr)
    tc_corr = np.interp(wet_bulb,x_wet_bulb,tc_wet_bulb)
    sc_corr  = griddata((x_sens_cool,y_sens_cool), z_sens_cool, (wet_bulb,dry_bulb), method='cubic',fill_value = -99)
    #
    # if the temperatures values are out of the table ( *) , go to the nearest correction factor
    #
#    if sc_corr  == -99:
#        sc_corr = griddata((x_sens_cool,y_sens_cool), z_sens_cool, (wet_bulb,dry_bulb), method='nearest')
    pow_corr = np.interp(wet_bulb,x_wet_bulb,power_wet_bulb)
    return tc_corr,1.0*sc_corr,pow_corr

#
#
def air_temp_heating_corr(dry_bulb,df):
    #
    # correction factors for air temperature cooling
    # dry_bulb  are in F
    #
    x_heat = df['EA_DBT'].values
    hc_air_temp = df['HC_CORR'].values
    pow_air_temp = df['POW_CORR'].values
    hc_corr = np.interp(dry_bulb,x_heat,hc_air_temp)
    pow_corr = np.interp(dry_bulb,x_heat,pow_air_temp)
    return hc_corr,pow_corr
#

def wb_tc_corr(wet_bulb):
    y = 0.000669*wet_bulb**2-0.073066*wet_bulb+2.907533
    return y

def tc_gpm_corr(x_gpm):
    # comes from WAHPcorrector ( Kavanaugh)
    y = 0.0281*(x_gpm)+0.9715
    return y

def hc_gpm_corr(x_gpm):
    y = -0.0962*(x_gpm)**2+0.2878*(x_gpm)+0.8085
    return y

def powc_gpm_corr(x_gpm) :
    y = 0.1102*(x_gpm)**2-0.3157*(x_gpm)+1.206
    return y

def powh_gpm_corr(x_gpm) :
    y  =-0.0154*(x_gpm)**2+0.067*(x_gpm)+0.9483
    return y

def SST(cfm,TC):
    y = (0.00978*(cfm/TC)+ 0.42087)
    return y

def WA_heat_pump(EWT_S,DEB_S,CFM_L,mode,df):

    if mode == 'heating':
        TC_heating = df['HCH'].values
        TC_heating = TC_heating[np.isfinite(TC_heating)]
        Power_heating = df['POWH'].values
        Power_heating = Power_heating[np.isfinite(Power_heating)]
        COP = df['COP'].values
        COP = COP[np.isfinite(COP)]
        EWTv = df['EWT_SH'].values
        x = np.unique(EWTv)
        x = x[np.isfinite(x)]
        EWTmin = min(x)
        EWTmax = max(x)
        DEBSv = df['FLOW_SH'].values
        y = np.unique(DEBSv)
        y = y[np.isfinite(y)]
        DEBSmin = min(y)
        DEBSmax = max(y)
        DEBLv = df['CFM_H'].values
        z = np.unique(DEBLv)
        z = z[np.isfinite(z)]
        DEBLmin = min(z)
        DEBLmax = max(z)
        if EWT_S < EWTmin or EWT_S > EWTmax:
            print('temperature source doit etre',  EWTmin,'  et ',EWTmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if DEB_S < DEBSmin or DEB_S > DEBSmax:
            print('Le debit source doit etre',  DEBSmin,'  et ',DEBSmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if CFM_L < DEBLmin or CFM_L > DEBLmax:
            print('Le debit charge doit etre',  DEBLmin,'  et ',DEBLmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        [X,Y,Z] = np.meshgrid(x,y,z)
        N1 = len(x)
        N2 = len(y)
        N3 = len(z)
        V1 = np.ones((N1,N2,N3))
        V2 = np.ones((N1,N2,N3))
        V3 = np.ones((N1,N2,N3))
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V1[i1,i2,i3] = TC_heating[k]
                    k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V2[i1,i2,i3] = Power_heating[k]
                    k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V3[i1,i2,i3] = COP[k]
                    k = k+1
        inp = (x,y,z)
        val = (EWT_S,DEB_S,CFM_L)
        tc  = interpn(inp,V1,val,'linear')[0]
        power = interpn(inp,V2,val,'linear')[0]
        cop = interpn(inp,V3,val,'linear')[0]
        return tc,power,cop
    else:
        TC_cooling = df['TCC'].values
        SC_cooling = df['SCC'].values
        Power_cooling = df['POWC'].values
        EER = df['EER'].values
        EWTv = df['EWT_SC'].values
        x = np.unique(EWTv)
        EWTmin = min(x)
        EWTmax = max(x)
        DEBSv = df['FLOW_SC'].values
        y = np.unique(DEBSv)
        DEBSmin = min(y)
        DEBSmax = max(y)
        DEBSl = df['CFM_C'].values
        z = np.unique(DEBSl)
        DEBLmin = min(z)
        DEBLmax = max(z)
        if EWT_S < EWTmin or EWT_S > EWTmax:
            print('temperature source doit etre',  EWTmin,'  et ',EWTmax)
            TC_cap =  -999
            SC_cap =  -999
            Pow =  -999
            EER =  -999
            return TC_cap,SC_cap,Pow,EER
        if DEB_S < DEBSmin or DEB_S > DEBSmax:
            print('Le debit source doit etre',  DEBSmin,'  et ',DEBSmax)
            TC_cap =  -999
            SC_cap =  -999
            Pow =  -999
            EER =  -999
            return TC_cap,SC_cap,Pow,EER
        if CFM_L < DEBLmin or CFM_L > DEBLmax:
            print('Le debit charge doit etre',  DEBLmin,'  et ',DEBLmax)
            TC_cap =  -999
            SC_cap =  -999
            Pow =  -999
            EER =  -999
            return TC_cap,SC_cap,Pow,EER
        N1 = len(x)
        N2 = len(y)
        N3 = len(z)
        V1 = np.ones((N1,N2,N3))
        V2 = np.ones((N1,N2,N3))
        V3 = np.ones((N1,N2,N3))
        V4 = np.ones((N1,N2,N3))
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V1[i1,i2,i3] = TC_cooling[k]
                    k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V2[i1,i2,i3] = SC_cooling[k]
                    k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V3[i1,i2,i3] = Power_cooling[k]
                    k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    V4[i1,i2,i3] = EER[k]
                    k = k+1
        inp = (x,y,z)
        val = (EWT_S,DEB_S,CFM_L)
        tc  = interpn(inp,V1,val,'linear')[0]
        sc  = interpn(inp,V2,val,'linear')[0]
        power = interpn(inp,V3,val,'linear')[0]
        eer = interpn(inp,V4,val,'linear')[0]
        return tc,sc,power,eer




def WW_heat_pump(EWT_S,DEB_S,EWT_L,DEB_L,mode,df):

    if mode == 'heating':
        TC_heating = df['HCH'].values
        TC_heating = TC_heating[np.isfinite(TC_heating)]
        Power_heating = df['POWH'].values
        Power_heating = Power_heating[np.isfinite(Power_heating)]
        COP = df['COP'].values
        COP = COP[np.isfinite(COP)]
        EWTSv = df['EWT_SH'].values
        x = np.unique(EWTSv)
        x = x[np.isfinite(x)]
        DEBSv = df['FLOW_SH'].values
        y = np.unique(DEBSv)
        y = y[np.isfinite(y)]
        EWTLv = df['EWT_LH'].values
        z = np.unique(EWTLv)
        z = z[np.isfinite(z)]
        DEBLv = df['FLOW_LH'].values
        w = np.unique(DEBLv)
        w = w[np.isfinite(w)]
        EWTSmin = min(x)
        EWTSmax = max(x)
        DEBSmin = min(y)
        DEBSmax = max(y)
        EWTLmin = min(z)
        EWTLmax = max(z)
        DEBLmin = min(w)
        DEBLmax = max(w)
        if EWT_S < EWTSmin or EWT_S > EWTSmax:
            print('temperature source doit etre',  EWTSmin,'  et ',EWTSmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if EWT_L < EWTLmin or EWT_L > EWTLmax:
            print('temperature source doit etre',  EWTLmin,'  et ',EWTLmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if DEB_S < DEBSmin or DEB_S > DEBSmax:
            print('Le debit source doit etre',  DEBSmin,'  et ',DEBSmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if DEB_L < DEBLmin or DEB_L > DEBLmax:
            print('Le debit charge doit etre',  DEBLmin,'  et ',DEBLmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        [X,Y,Z,W] = np.meshgrid(x,y,z,w)
        N1 = len(x)
        N2 = len(y)
        N3 = len(z)
        N4 = len(w)
        V1 = np.ones((N1,N2,N3,N4))
        V2 = np.ones((N1,N2,N3,N4))
        V3 = np.ones((N1,N2,N3,N4))
        V4 = np.ones((N1,N2,N3,N4))
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V1[i1,i2,i3,i4] = TC_heating[k]
                        k = k+1
        k = 0
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V3[i1,i2,i3,i4] = Power_heating[k]
                        k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V4[i1,i2,i3,i4] = COP[k]
                        k = k+1
        inp = (x,y,z,w)
        val = (EWT_S,DEB_S,EWT_L,DEB_L)
        TC_cap  = interpn(inp,V1,val,'linear')[0]
        Pow = interpn(inp,V3,val,'linear')[0]
        cop = interpn(inp,V4,val,'linear')[0]
        return TC_cap,Pow,cop
    else:
        TC_cooling = df['TC'].values
        TC_cooling = TC_cooling[np.isfinite(TC_cooling)]
        Power_cooling = df['POWC'].values
        Power_cooling = Power_cooling[np.isfinite(Power_cooling)]
        EER = df['EER'].values
        EER = EER[np.isfinite(EER)]
        EWTSv = df['EWT_SC'].values
        x = np.unique(EWTSv)
        x = x[np.isfinite(x)]
        DEBSv = df['FLOW_SC'].values
        y = np.unique(DEBSv)
        y = y[np.isfinite(y)]
        EWTLv = df['EWT_LC'].values
        z = np.unique(EWTLv)
        z = z[np.isfinite(z)]
        DEBLv = df['FLOW_LC'].values
        w = np.unique(DEBLv)
        w = w[np.isfinite(w)]
        EWTSmin = min(x)
        EWTSmax = max(x)
        DEBSmin = min(y)
        DEBSmax = max(y)
        EWTLmin = min(z)
        EWTLmax = max(z)
        DEBLmin = min(w)
        DEBLmax = max(w)
        if EWT_S < EWTSmin or EWT_S > EWTSmax:
            print('temperature source doit etre',  EWTmin,'  et ',EWTmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if EWT_L < EWTLmin or EWT_L > EWTLmax:
            print('temperature source doit etre',  EWTmin,'  et ',EWTmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if DEB_S < DEBSmin or DEB_S > DEBSmax:
            print('Le debit source doit etre',  DEBSmin,'  et ',DEBSmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        if DEB_L < DEBLmin or DEB_L > DEBLmax:
            print('Le debit charge doit etre',  DEBLmin,'  et ',DEBLmax)
            TC_cap =  -999
            Pow =  -999
            COP =  -999
            return TC_cap,Pow,COP
        [X,Y,Z,W] = np.meshgrid(x,y,z,w)
        N1 = len(x)
        N2 = len(y)
        N3 = len(z)
        N4 = len(w)
        V1 = np.ones((N1,N2,N3,N4))
        V2 = np.ones((N1,N2,N3,N4))
        V3 = np.ones((N1,N2,N3,N4))
        V4 = np.ones((N1,N2,N3,N4))
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V1[i1,i2,i3,i4] = TC_cooling[k]
                        k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V3[i1,i2,i3,i4] = Power_cooling[k]
                        k = k+1
        k = 0
        for i1 in range(0,N1):
            for i2 in range(0,N2):
                for i3 in range(0,N3):
                    for i4 in range(0,N4):
                        V4[i1,i2,i3,i4] = EER[k]
                        k = k+1
        inp = (x,y,z,w)
        val = (EWT_S,DEB_S,EWT_L,DEB_L)
        TC_cap  = interpn(inp,V1,val,'linear')[0]
        Pow = interpn(inp,V3,val,'linear')[0]
        eer = interpn(inp,V4,val,'linear')[0]
        return TC_cap,Pow,eer


def WA_hp_psi(EWT_S,DEB_S,df):
    try:
        PSIv = df['PSI'].values
    except Exception as e:
        type(e).__name__ + ': ' + str(e)
        erru = type(e).__name__
        print('There is no column named '+ str(e) + 'in the file')
        return -999
    PSIv = PSIv[np.isfinite(PSIv)]
    EWTv = df['EWT_SC'].values
    DEBSv = df['FLOW_SH'].values
    DEBLv = df['CFM_H'].values
    EWTv = EWTv[np.isfinite(EWTv)]
    x = np.unique(EWTv)
    nt = len(x)
    x = x[np.isfinite(x)]
    EWTmin = min(x)
    EWTmax = max(x)
    DEBSv = DEBSv[np.isfinite(DEBSv)]
    y = np.unique(DEBSv)
    y = y[np.isfinite(y)]
    z = np.unique(DEBLv)
    z = z[np.isfinite(z)]
    nx = len(z)
    nd = len(y)
    nz = nd*nx
    it = cherche_index(EWT_S,x)
    if it < nt-1:
        i1 = it
        i2 = it+1
        i3 = it+2
        z1 = PSIv[i1*nz:i2*nz:nx]
        z2 = PSIv[i2*nz:i3*nz:nx]
        zn = np.zeros(nd)
        for i in range(0,nd):
            zn[i] = z1[i] + (z2[i] - z1[i])/(x[i2]-x[i1])*(EWT_S - x[i1])
        p = np.polyfit(y,zn,2)
        psi = np.polyval(p,DEB_S)
        return psi
    else:
        i1 = it
        i2 = it+1
        zn = PSIv[i1*nz:i2*nz:nx]
        p = np.polyfit(y,zn,2)
        psi = np.polyval(p,DEB_S)
        return psi

def HP_flow_rates(df,mode = 'heating'):
    if mode == 'heating':
        DEBSv = df['FLOW_SH'].values
    else:
        DEBSv = df['FLOW_SC'].values
    DEBSv = DEBSv[np.isfinite(DEBSv)]
    x = np.unique(DEBSv)
    gpmm = np.median(x)
    gpm_min = np.min(x)
    gpm_max = np.max(x)
    return gpmm,gpm_min,gpm_max

def HP_cfm_rates(df,mode = 'heating'):
    if mode == 'heating':
        DEBSv = df['CFM_H'].values
    else:
        DEBSv = df['CFM_C'].values
    DEBSv = DEBSv[np.isfinite(DEBSv)]
    x = np.unique(DEBSv)
    cfm_min = np.min(x)
    cfm_max = np.max(x)
    return cfm_min,cfm_max