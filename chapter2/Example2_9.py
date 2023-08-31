import numpy as np
from scipy.optimize import newton,bisect
# exemple 2.9
A = 12.0
Ti = 22.0
To = -10.0
Tenv1 = -10.0
Tenv2 = -20.0
DT = Ti-To
ep = 0.9
sig = 5.67e-8
ho = 25.0
hi = 8.0
L1 = 0.24
k1 = 2.1
L2 = 0.16
k2 = 0.04
R1 = 1.0/(A*ho)
R2 = L1/(A*k1)
R3 = L2/(A*k2)
R4 = 1.0/(A*hi)
Rpp1 = 1.0/ho
Rpp2 = L1/k1
Rpp3 = L2/k2
Rpp4 = 1.0/hi
Tski = -9.7+273.15 # hypothese
Tok = To + 273.15
Tik = Ti + 273.15
Rppint = Rpp2+Rpp3+Rpp4
#
# a)
#
qsol = 0
Tenvk = Tenv1 + 273.15
hrad = ep*sig*(Tski+Tenvk)*(Tski**2+Tenvk**2)
Rpprad = 1.0/hrad
Rppext = 1.0/(ho+hrad)
Rpptot = Rppint + Rppext
Rtot = Rpptot/A
qa = DT/Rtot
qppa = qa/A
U = 1.0/Rpptot
Tsa = Ti - qppa*Rppint
qppco = (Tsa-To)*ho
qppra = (Tsa-Tenv1)*hrad
qa = qppa*A
U = qppa/DT
print('Solution using the approxiamte radiative resistance')
print ('a) Heat losses =  ',qppa, ' W/m2')
print ('a) Heat losses =   ',qa, ' W')
print ('a) Global loss coefficient U  =   ',U ,' W/m2 K')
#
# exact solution
#
def fct(Tsk):
    qout = (Tsk-Tok)/Rpp1 + ep*sig*(Tsk**4-Tenvk**4)
    qin = qsol + (Tik - Tsk)/Rppint
    y = qout - qin
    return y
#
Tsn = 250 # hypothese pour le processus itératif
#
#
Tssk1  = newton(fct,Tsn)
Tssk = Tssk1
print ('Tsa (new ) = ',Tssk - 273.15)
hrada = ep*sig*(Tssk+Tenvk)*(Tssk**2+Tenvk**2)
Rpptot2 = Rppint + 1.0/(ho+hrada)
Rtot2 = Rpptot2/A
Tsa2 = Tssk - 273.15
qppa2 = (Ti-Tsa2)/Rppint
qppaco2 = (Tsa2-To)/Rpp1
qppara2 = ep*sig*(Tssk**4-Tenvk**4)
qa2 = qppa2*A
U2 = qppa2/DT
print('Solution using the exact expression')
print ('a) Heat losses =  ',qppa2, ' W/m2')
print ('a) Heat losses =   ',qa2, ' W')
print ('a) Global loss coefficient U  =   ',U2 ,' W/m2 K')


#
# b)
# solution approchée
#
Tenvk = Tenv2 + 273.15
hradb = ep*sig*(Tski+Tenvk)*(Tski**2+Tenvk**2)
Rpprad = 1 /hradb
Tsb = (To/Rpp1 + Ti/Rppint + Tenv2/Rpprad)/(1.0/Rpp1 + 1.0/Rppint + 1.0/Rpprad)
Tsbk = Tsb + 273.15
hradbn = ep*sig*(Tski+Tenvk)*(Tski**2+Tenvk**2)
qppb = (Ti-Tsb)/Rppint
qppbco = (Tsb-To)/Rpp1
qppbra = (Tsb-Tenv2)/Rpprad
qppo = qppbco + qppbra
qb = qppb*A
Ub = qppb/DT
print('Solution using the approxiamte radiative resistance')
print ('Tsb  = ',Tsbk - 273.15)
print ('b) Heat losses =  ',qppb, ' W/m2')
print ('b) Heat losses =   ',qb, ' W')
print ('b) Global loss coefficient U  =   ',Ub ,' W/m2 K')

#
# solution exacte
#
Tenvk = Tenv2 + 273.15
Tsskb  = newton(fct,Tski)
hradb = ep*sig*(Tsskb+Tenvk)*(Tsskb**2+Tenvk**2)
Tsb2 = Tsskb - 273.15
qppb2 = (Ti-Tsb2)/Rppint
qppbco2 = (Tsb2-To)/Rpp1
qppbra2 = ep*sig*(Tsskb**4-Tenvk**4)
qb2 = qppb2*A
Tssb = Ti - qppb2*Rppint
Ub2 = qppb2/DT
print('Solution using the exact expression')
print ('Tsb  = ',Tsskb - 273.15)
print ('b) Heat losses =  ',qppb2, ' W/m2')
print ('b) Heat losses =   ',qb2, ' W')
print ('b) Global loss coefficient U  =   ',Ub2 ,' W/m2 K')

#
# c)
#
# solution
#
Tenvk = Tenv1 + 273.15
qsol = 0.8*750
ok = False
compt = 0
Tsk = Tski            # hypothese initial Ts = -9.7 ( Exemple 2.1 )
Rppext = 1.0/(ho+hrad)
Tsc = (To/Rppext + Ti/Rppint + qsol)/(1.0/Rppext + 1.0/Rppint)
qppci = (Ti-Tsc)/Rppint
print ('Tsc  = ',Tsc)
print ('c) Heat losses =  ',qppci, ' W/m2')

while not ok:
    hradc = ep*sig*(Tsk+Tenvk)*(Tsk**2+Tenvk**2)
#    print('hradc = ',hradc)
    Rpprad = 1.0/hradc
    Rppext = 1.0/(ho+hradc)
    Tsc = (To/Rppext + Ti/Rppint + qsol)/(1.0/Rppext + 1.0/Rppint)
    qppc = (Ti-Tsc)/Rppint
    qppext = (Tsc-To)/Rppext
    qc = qppc*A
    Tskn = Tsc + 273.15
    if abs(Tskn - Tsk) < 0.001:
        ok = True
    else:
        Tsk = Tskn
        compt = compt + 1
print ('Tsc  = ',Tskn)
print ('Tsc  = ',Tskn - 273.15)
print ('c) Heat losses =  ',qppc, ' W/m2')
print ('c) Heat losses =   ',qc, ' W')

#
# solution exacte
#
Tsskc  = newton(fct,Tsk)
Tsc2 = Tsskc - 273.15
print ('Tsc  = ',Tsskc)
print ('Tsc  = ',Tsc2)
qppc2 = (Ti-Tsc2)/Rppint
qppcco2 = (Tsc2-To)/Rpp1
qppcra2 = ep*sig*(Tsskc**4-Tenvk**4)
qppextc2 =qppcco2 + qppcra2
qc2 = qppc2*A
print ('c) Heat losses =  ',qppc2, ' W/m2')
print ('c) Heat losses =   ',qc2, ' W')

# d
#
# solution
#
Tenvk = Tenv2 + 273.15
qsol = 0.8*750
ok = False
compt = 0
Tsk = Tski            # hypothese initial Ts = -9.7 ( Exemple 2.1 )
Tsda = (To/Rpp1 + Ti/Rppint + Tenv2/Rpprad + qsol)/(1.0/Rpp1 + 1.0/Rppint + 1.0/Rpprad)
qppda = (Ti-Tsda)/Rppint

Tssk1  = newton(fct,Tsn)
Tsdb = Tssk1 - 273.15
qppdb = (Ti-Tsdb)/Rppint
