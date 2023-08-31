#coding: latin-1
# exemple 3.3
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
import pandas as pd
from conversion_md import *
from geothermal_md import *
from heat_pump_md import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
from scipy.interpolate import interp1d
#
# data
#
patm = 101.325*1000.0
Tdbc = 25.0         # air dry bulb temperaure
Tdbk = Tdbc + 273.15
phi = 0.4          # air humidity
ewt = 30.0          # entering water temperaure
cfm = 1200        # cfm of air in air side HP
cfm_nominal = 1300.0
gpm_nom = 8
gpm = 6

#
#wet bulb temperaure
#
Twbk = HAPropsSI('Twb','T',Tdbk,'P',patm,'RH',phi)
Twbc = Twbk  - 273.15
print ('Wet-bulb = ',Twbc )
EWTF = F_C(ewt)
dry_bulb = F_C(Tdbc)
wet_bulb = F_C(Twbc)
heat_pump = '..\\data\\wa_nv036.xls'
df = pd.read_excel(heat_pump)
gpmv,gpm_min,gpm_max = HP_flow_rates(df)
cfm_loadv  = df['CFM_C'].values
cfm_loadv = cfm_loadv[np.isfinite(cfm_loadv)]
cfm_min = min(cfm_loadv)
cfm_max = max(cfm_loadv)
TC,SC,Pow,EER = WA_heat_pump(EWTF,gpm,cfm,'cooling',df)
print ('nominal total capacity = ',TC, ' kbtu/hr' )
print ('nominal power = ',Pow, ' kW')
print ('nominal EER = ',EER)
hp_comp = 'wf'
file2 = '..\\data\\'+ hp_comp +'_cc_temp_corr.xls'
df_temp = pd.read_excel(file2)
F_tc_aT,F_sc_aT,F_pow_aT = air_temp_cooling_corr(dry_bulb,wet_bulb,df_temp)
TC2 = TC*F_tc_aT
SC2 = SC*F_sc_aT
Pow2 = Pow*F_pow_aT
print ('corrected total capacity = ',TC2, ' kbtu/hr' )
print ('corrected sensible capacity = ',SC2, ' kbtu/hr' )
print ('corrected power = ',Pow2, ' kW')
print ('corrected eer = ',TC2/Pow2)
