#coding: latin-1
import numpy as np

def shell_tube_eff(NTU,Cr,nshell=1):
#
# shell and tube (nshell shell pass)
# NTU is the total NTU
    NTUn = NTU/nshell
    G = np.sqrt(1+Cr**2)
    y = np.exp(-NTUn*G)
    ep1 = 2/(1+Cr+G*(1+y)/(1-y))
    if nshell > 1:
        if Cr == 1:
            ep = nshell*ep1/(1+ep1*(n-1))
        else:
            z = (1-ep1*Cr)/(1-ep1)
            ep = (z**nshell-1)/(z**nshell-Cr)
    else:
        ep = ep1
    return ep

def shell_tube_NTU(ep,Cr,nshell=1):
#
# shell and tube (nshell shell pass)
# NTU is the total NTU
    G = np.sqrt(1+Cr**2)
    if nshell > 1:
        if Cr ==1:
            ep1 = ep/(nshell - ep*(n-1))
        else:
            F = ((ep*Cr -1)/(ep-1))**(1/nshell)
            ep1 = (F-1)/(F-Cr)
    else:
        ep1 = ep
    E = (2/ep1 - (1+Cr))/G
    if E > 1:
        NTU1 = -np.log((E-1)/(E+1))/G
        NTU = nshell*NTU1
    else:
        print('impossible')
        NTU = -999
    return NTU
def F_coef_shell_tube(Tci,Tco,Thi,Tho,N):
    P = (Tco-Tci)/(Thi-Tci)
    R = (Thi-Tho)/(Tco-Tci)
    A = 2/P - 1 - R
    B = 2/P*np.sqrt((1-P)*(1-P*R))
    num = np.sqrt(R**2 + 1)/(R-1)*np.log((1-P)/(1-P*R))
    if N == 1:
        den = np.log((A+np.sqrt(R**2+1))/(A-np.sqrt(R**2+1)))
        F = num/den
    else:
        den = 2*np.log((A+B+np.sqrt(R**2+1))/(A+B-np.sqrt(R**2+1)))
        F = num/den
    return F


def counter_flow_NTU(ep,Cr):
    if Cr < 1:
        NTU = 1/(Cr-1)*np.log((ep-1)/(ep*Cr-1))
    else:
        NTU = ep/(1.0-ep)
    return NTU

def counter_flow_eff(NTU,Cr):

    if Cr < 1:
        ep = (1-np.exp(-NTU*(1-Cr)))/(1-Cr*np.exp(-NTU*(1-Cr)))
    else:
        ep = NTU/(1+NTU)
    return ep