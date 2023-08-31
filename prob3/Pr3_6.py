import numpy as np


def calcul_CAP(T_L,T_S):
    y =  66 -0.52*T_L + 2.74*T_S  -0.019*T_L*T_S - 0.023*T_S**2
    return y
def calcul_COP(T_L,T_S):
    y =  3.78 - 0.042*T_L + 0.076*T_S  -0.00086*T_L*T_S - 0.00018*T_S**2
    return y
mph = 200./60.
Cph = 4180
Cpc = 4180
Ch = mph*Cph
mpc = 70.0/60.0
ef_ech   = 0.8

# 1st case
T_load = 10
T_source =  30
T_hp_o = 65
T_hp_in = T_load
T_ec_in = T_load
cap1 = calcul_CAP(T_hp_in,T_source)
cop1 = calcul_COP(T_hp_in,T_source)
cap1 = cap1*1000   # W
mp_hp1 = cap1/(Cpc*(T_hp_o-T_hp_in))
pow1 = cap1/cop1
mp_ec1 = mpc - mp_hp1
T_source_oa1 = T_source - cap1/Ch
#print('T source out hp',T_source_oa1)
# heat exchanger
Cc = mp_ec1*Cpc
Cmin =  min(Cc,Ch)
qmax = Cmin*(T_source_oa1 - T_ec_in)
qech1 = ef_ech*qmax
T_source_ob1 = T_source_oa1 - qech1/Ch
#print('Tsource out ech = ', T_source_ob1)
T_ec_o1 = T_ec_in + qech1/Cc
T_load_o1a = (mp_ec1*T_ec_o1 + mp_hp1*T_hp_o)/mpc
q_tot1 = qech1 + cap1
copf1 = q_tot1/pow1
T_load_o1b = T_load + q_tot1/(mpc*Cpc)
#print('Tout = ', T_load_o1a)
#print('Tout = ', T_load_o1b)
print('qtot 1 = ', q_tot1/1000)
print('cop sys 1 = ', copf1)

# 2st case
ok = False
compt = 1
compt_max = 100
delt = 0.001
T_s_in = T_source_oa1
while not ok:
    cap2 = calcul_CAP(T_hp_in,T_s_in)
    cap2 = cap2*1000   # W
    mp_hp2 = cap2/(Cpc*(T_hp_o-T_hp_in))
    mp_ec2 = mpc - mp_hp2
    Cc = mp_ec2*Cpc
    Cmin =  min(Cc,Ch)
    qmax = Cmin*(T_source - T_ec_in)
    qech2 = ef_ech*qmax
    T_source_oa2 = T_source - qech2/Ch
    dif = abs(T_source_oa2 -T_s_in)
    if dif > delt:
        T_s_in =  T_source_oa2
        compt = compt + 1
        if compt > compt_max:
            print('erreur')
            ok = True
    else:
        ok = True
cop2 = calcul_COP(T_hp_in,T_source_oa2)
pow2 = cap2/cop2
T_source_ob2 = T_source_oa2 - cap2/Ch
#print('T source out ech',T_source_oa2)
#print('T source out hp',T_source_ob2)
T_ec_o2 = T_ec_in + qech2/Cc
T_load_o2a = (mp_ec2*T_ec_o2 + mp_hp2*T_hp_o)/mpc
q_tot2 = qech2 + cap2
copf2 = q_tot2/pow2
T_load_o2b = T_load + q_tot2/(mpc*Cpc)
#print('Tout = ', T_load_o2a)
#print('Tout = ', T_load_o2b)
print('qtot 2 = ', q_tot2/1000)
print('cop sys 2 = ', copf2)

#
# 3st case
#
mp_ec3 = mpc
Cc = mp_ec3*Cpc
Cmin =  min(Cc,Ch)
qmax = Cmin*(T_source - T_ec_in)
qech3 = ef_ech*qmax
T_source_oa3 = T_source - qech3/Ch
T_ec_o3 = T_ec_in + qech3/Cc
T_hp_in3 = T_ec_o3
cap3 = calcul_CAP(T_hp_in3,T_source_oa3)
cap3 = cap3*1000   # W
mp_hp3 = cap3/(Cpc*(T_hp_o-T_hp_in3))
cop3 = calcul_COP(T_hp_in3,T_source_oa3)
pow3 = cap3/cop3
T_source_ob3 = T_source_oa3 - cap3/Ch
#print('T source out ech',T_source_oa3)
#print('T source out hp',T_source_ob3)
mp_bp = mpc - mp_hp3
T_load_o3a = (mp_bp*T_ec_o3 + mp_hp3*T_hp_o)/mpc
q_tot3 = qech3 + cap3
copf3 = q_tot3/pow3
T_load_o3b = T_load + q_tot3/(mpc*Cpc)
#print('Tout = ', T_load_o3a)
#print('Tout = ', T_load_o3b)
print('qtot 3 = ', q_tot3/1000)
print('cop sys  3 = ', copf3)

#
# 4st case
#
mp_ec4 = mpc
Cc = mp_ec4*Cpc
Cmin =  min(Cc,Ch)
ok = False
compt = 1
compt_max = 100
delt = 0.001
T_hp_in4 = T_ec_o3
T_ec_in4 = T_load
while not ok:
    cap4 = calcul_CAP(T_hp_in4,T_source)
    cap4 = cap4*1000   # W
    mp_hp4 = cap4/(Cpc*(T_hp_o-T_hp_in4))
    T_source_oa4 = T_source - cap4/Ch
    qmax = Cmin*(T_source_oa4 - T_ec_in4)
    qech4 = ef_ech*qmax
    T_ec_o4 = T_ec_in4 + qech4/Cc
    dif = abs(T_ec_o4 -T_hp_in4)
    if dif > delt:
        T_hp_in4 =  T_ec_o4
        compt = compt + 1
        if compt > compt_max:
            print('erreur')
            ok = True
    else:
        ok = True
cop4 = calcul_COP(T_hp_in4,T_source)
pow4 = cap4/cop4
T_source_ob4 = T_source_oa4 - qech4/Ch
#print('T source out hp',T_source_oa4)
#print('T source out ech',T_source_ob4)
mp_bp4 = mpc - mp_hp4
T_load_o4a = (mp_bp4*T_ec_o4 + mp_hp4*T_hp_o)/mpc
q_tot4 = qech4 + cap4
copf4 = q_tot4/pow4
T_load_o4b = T_load + q_tot4/(mpc*Cpc)
#print('Tout = ', T_load_o4a)
#print('Tout = ', T_load_o4b)
print('qtot 4= ', q_tot4/1000)
print('cop sys 4= ', copf4)


