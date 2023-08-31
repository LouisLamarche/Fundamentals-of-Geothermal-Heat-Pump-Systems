#coding: utf-8
#-------------------------------------------------------------------------------
from numpy import *
from  matplotlib.pyplot import *
import pandas as pd
donnees = loadtxt("..\data\loads_building1.txt")   # charges horaires en 2 colonnes, chauffage et climatisation pour chaque heure de l'ann?e
qch_hr = donnees[:,0]   # charges chauffage en kW
qcl_hr = donnees[:,1]   # charges climatisation en kW
qtot = qch_hr-qcl_hr    # charge horaires totales + chauffage , - clim
hrma = array([744,672,744,720,744,720,744,744,720,744,720,744])
hrmb = 730*ones(12,dtype=int)
Qm_cha = ones(12)    # consommations  mensuelles chauffage
Qm_cla = ones(12)    # consommations  mensuelles climatisation
Qm_chb = ones(12)    # consommations  mensuelles chauffage
Qm_clb = ones(12)    # consommations  mensuelles climatisation
Qm = ones(12)    # consommations  mensuelles climatisation
p_ch = ones(12)     # puissance maximales mensuelles chauffage
p_cl = ones(12)     # puissance maximales mensuelles climatisation
pm_ch = ones(24)     # puissance moyenne mensuelles chauffage
pm_cl = ones(24)     # puissance moyennes mensuelles climatisation
# vecteur pour la représentation graphique
pu_ch = ones(24)
pu_cl = ones(24)
Qu_ch = ones(24)
Qu_cl = ones(24)
i1 = 0
for i in range (0,12):
    i2 = i1+hrma[i]-1
    Qm_cha[i] = sum(qch_hr[i1:i2+1])
    Qm_cla[i] = sum(qcl_hr[i1:i2+1])
    Qm[i] = sum(qtot[i1:i2+1])
    p_ch[i] = mean(qch_hr[i1:i2+1])
    p_cl[i] = mean(qcl_hr[i1:i2+1])
    pu_cl[2*i:2*i+2] = -p_cl[i]
    Qu_cl[2*i:2*i+2] = Qm_cla[i]
    pu_ch[2*i:2*i+2] = p_ch[i]
    Qu_ch[2*i:2*i+2] = Qm_cha[i]
    pm_cl[2*i:2*i+2] = Qm_cla[i]/hrma[i]
    pm_ch[2*i:2*i+2] = Qm_cha[i]/hrma[i]
    i1 = i2+1
i1 = 0
for i in range (0,12):
    i2 = i1+hrmb[i]-1
    Qm_chb[i] = sum(qch_hr[i1:i2+1])
    Qm_clb[i] = sum(qcl_hr[i1:i2+1])
    i1 = i2+1

t = arange(0,8760)
im = arange(0,12)
mm = zeros(24)
for i in range(0,11):
    mm[2*i+1:2*i+3] = i+1
mm[23] = 12

# représentation graphique
#
flag_plot = 0
if flag_plot == 1:
    fig1 = figure()
    plot(mm,Qu_ch,linewidth=2)
    ylabel('kWh',fontsize = 14)
    xlabel('month',fontsize = 14)
    title('Heating',fontsize = 16)
    fig2 = figure()
    plot(mm,Qu_cl,linewidth=2)
    ylabel('kWh',fontsize = 14)
    xlabel('month',fontsize = 14)
    title('Cooling',fontsize = 16)
    show()
else:
    data =  np.zeros((12,4))
    for ip in range(0,12):
        data[ip,0] = Qm_cha[ip]
        data[ip,1] = Qm_cla[ip]
        data[ip,2] = Qm_chb[ip]
        data[ip,3] = Qm_clb[ip]
    indd  = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    coll = ['Heating (kWh)','Cooling (kWh)','Heating (kWh)','Cooling (kWh)']
    df = pd.DataFrame(data,columns = coll,index = indd)
    def f1(x):
        return '%.0f' % x
    print(df.to_latex(index=True,formatters=[f1,f1,f1,f1]))


