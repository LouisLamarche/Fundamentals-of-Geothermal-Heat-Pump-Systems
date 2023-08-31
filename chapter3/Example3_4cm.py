#coding: utf-8
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
import pandas as pd
from conversion_md import *
from geothermal_md import *
from heat_pump_md import *
from CoolProp.CoolProp import *
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
cfm_nominal = 1000.0
gpm_nom = 5.6
gpm = 6
#
# heat pump
heat_pump = 'cm'   # climate master tranquility
#heat_pump = 'wf'   # Water furnace series 7

file1 = '..\\data\\' + heat_pump+'_air_flow_corr.xls'
df_flow = pd.read_excel(file1)
file2 = '..\\data\\' + heat_pump +'_cc_temp_corr.xls'
df_temp = pd.read_excel(file2)
file3 = '..\\data\\' + heat_pump +'_hc_temp_corr.xls'
df_hc = pd.read_excel(file3)

#
#wet bulb temperaure
#
Twbk = HAPropsSI('Twb','T',Tdbk,'P',patm,'RH',phi)
Twbc = Twbk  - 273.15
EWTF = F_C(ewt)
Tdbf = F_C(Tdbc)
print ('Dry-bulb = ',Tdbf,' F')
Twbf = F_C(Twbc)
print ('Wet-bulb = ',Twbf,' F')
#
# correction for water gpm
#
x_gpm = gpm/gpm_nom
F_tc_gpm = tc_gpm_corr(x_gpm)
F_powc = powc_gpm_corr(x_gpm)
F_hc = hc_gpm_corr(x_gpm)
F_powh = powh_gpm_corr(x_gpm)
xt = np.array([59,77,86])
tct = np.array([30.8,29.2,28])
eert = np.array([26.7,19.4,17.3])
f_tc = interp1d(xt,tct,'linear',fill_value = 'extrapolate')
f_eer = interp1d(xt,eert,'linear',fill_value = 'extrapolate')
tcn = f_tc(EWTF)
tc77  = f_tc(77)
eern = f_eer(EWTF)
pown = tcn/eern
print ('nominal total capacity = ',tcn, ' kbtu/hr' )
print ('nominal power = ',pown, ' kW')
#
# air flow correction
#
x_air = cfm/cfm_nominal
F_tc_cfm,F_sc_cfm,F_pow_cfm = air_flow_corr(x_air,df_flow)
dry_bulb = F_C(Tdbc)
wet_bulb = F_C(Twbc)
#
# air temp corr
#
F_tc_aT,F_sc_aT,F_pow_aT = air_temp_cooling_corr(dry_bulb,wet_bulb,df_temp)
# if F_sc_aT = -99,  SC = TC
tc1 = tcn*F_tc_gpm*F_tc_cfm*F_tc_aT
pow1 = pown*F_powc*F_pow_cfm*F_pow_aT
print ('corrected total capacity = ',tc1, ' kbtu/hr' )
print ('corrected power = ',pow1, ' kW')
print ('corrected eer = ',tc1/pow1)
