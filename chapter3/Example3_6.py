#coding: utf-8
#
# exemple 3.5
# written by Louis Lamarche
# 15 october 2017
#
import numpy as np
E = 25000.0  # consumption (KWh)
eff = 0.8
Ec = E/eff  # gaz consumption
GWP_oil = 0.295       # kg/KWh
GWP_gaz = 0.220       # kg/KWh
GWP_elec_hydro = 0.002  # kg/kWh
GWP_elec_oil = 0.960  # kg/kWh
GWP_elec_gaz = 0.58  # kg/kWh
GWP_elec_coal = 1.14  # kg/kWh
GES_oil  = GWP_oil*Ec
GES_gaz  = GWP_gaz*Ec
GES_elec_hydro  = GWP_elec_hydro*E
GES_elec_oil  = GWP_elec_oil*E
GES_elec_gaz  = GWP_elec_gaz*E
GES_elec_coal  = GWP_elec_coal*E
# Solution HP
# direct contribution
Capacite1 = 4.0      # kg
Capacite2 = 6.5    # kg
leak = 0.03   # hypothesis
GWP_410 = 1730.0
GES_direct1 = Capacite1*leak*GWP_410
GES_direct2 = Capacite2*leak*GWP_410
COP1 = 2.5
COP2 = 3.5
E1 = E/COP1
E2 = E/COP2
GES_elec_hydro1  = GWP_elec_hydro*E1
GES_elec_hydro2  = GWP_elec_hydro*E2
GES_elec_oil1  = GWP_elec_oil*E1
GES_elec_gaz1  = GWP_elec_gaz*E1
GES_elec_coal1  = GWP_elec_coal*E1
GES_elec_oil2  = GWP_elec_oil*E2
GES_elec_gaz2  = GWP_elec_gaz*E2
GES_elec_coal2  = GWP_elec_coal*E2
GES_PAC_AA_hydro = GES_direct1 + GES_elec_hydro1
GES_PAC_AA_oil = GES_direct1 + GES_elec_oil1
GES_PAC_AA_gaz = GES_direct1 + GES_elec_gaz1
GES_PAC_AA_coal = GES_direct1 + GES_elec_coal1
GES_PAC_GE_hydro = GES_direct1 + GES_elec_hydro2
GES_PAC_GE_oil = GES_direct1 + GES_elec_oil2
GES_PAC_GE_gaz = GES_direct1 + GES_elec_gaz2
GES_PAC_GE_coal = GES_direct1 + GES_elec_coal2
GES_PAC_DX_hydro = GES_direct2 + GES_elec_hydro2
GES_PAC_DX_oil = GES_direct2 + GES_elec_oil2
GES_PAC_DX_gaz = GES_direct2 + GES_elec_gaz2
GES_PAC_DX_coal = GES_direct2 + GES_elec_coal2

print ('Hydro')
print ('\t Gaz furnace ' + str(GES_gaz))
print ('\t Oil furnace ' + str(GES_oil))
print ('\t Electric heating  ' + str(GES_elec_hydro))
print ('\t Air-Air HP ' + str(GES_PAC_AA_hydro))
print ('\t SL-GSHP ' + str(GES_PAC_GE_hydro))
print ('\t DX-GSHP  ' + str(GES_PAC_DX_hydro))
print ('Coal')
print ('\t Gaz furnace ' + str(GES_gaz))
print ('\t Oil furnace ' + str(GES_oil))
print ('\t Electric heating  ' + str(GES_elec_coal))
print ('\t Air-Air HP ' + str(GES_PAC_AA_coal))
print ('\t SL-GSHP ' + str(GES_PAC_GE_coal))
print ('\t DX-GSHP  ' + str(GES_PAC_DX_coal))
print (' Gaz')
print ('\t Gaz furnace ' + str(GES_gaz))
print ('\t Oil furnace ' + str(GES_oil))
print ('\t Electric heating  ' + str(GES_elec_gaz))
print ('\t Air-Air HP ' + str(GES_PAC_AA_gaz))
print ('\t SL-GSHP ' + str(GES_PAC_GE_gaz))
print ('\t DX-GSHP  ' + str(GES_PAC_DX_gaz))
print (' oil')
print ('\t Gaz furnace ' + str(GES_gaz))
print ('\t Oil furnace ' + str(GES_oil))
print ('\t Electric heating  ' + str(GES_elec_oil))
print ('\t Air-Air HP ' + str(GES_PAC_AA_oil))
print ('\t SL-GSHP ' + str(GES_PAC_GE_oil))
print ('\t DX-GSHP  ' + str(GES_PAC_DX_oil))
#
#  Effect of selling electricity
#
Econ = E - E2
GESp = GES_gaz - GES_elec_hydro
GESm = GWP_elec_gaz*E
GES_sauve1 = GESp - GESm
print ('GAZ instead of baseboards ' + str(GES_sauve1))
GESp = GES_PAC_GE_hydro - GES_elec_hydro
GESm = GWP_elec_gaz*(E - E2)
GES_sauve2 = GESp - GESm
print ('GSHP instead of baseboards ' + str(GES_sauve2))
leak_tot1 = Capacite1*GWP_410
leak_tot2 = Capacite2*GWP_410
