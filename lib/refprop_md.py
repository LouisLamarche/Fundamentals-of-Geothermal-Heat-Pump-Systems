import numpy as np

from CoolProp.CoolProp import  *
from scipy.optimize import newton,brentq

set_config_string(ALTERNATIVE_REFPROP_PATH, 'c:\\Program Files (x86)\\REFPROP')
REFPROP = AbstractState('REFPROP','WATER&AMMONIA')
Mwater = PropsSI('WATER','molemass')
Mammonia = PropsSI('AMMONIA','molemass')

def calcul_state(T,p,x):
    REFPROP.set_mass_fractions([1 - x,x])
    REFPROP.update(PQ_INPUTS, p, 0)
    Tb = REFPROP.T()
    REFPROP.update(PQ_INPUTS, p, 1)
    Td = REFPROP.T()
    if T < Tb:
        etat = 'liquid'
    elif T > Td:
        etat = 'vapor'
    else:
        etat = 'mix'
    return etat



def change_ref(hi,x):
    dhrefa = 3.4787e+02  # kJ/kg
    dhrefw = -25.4596   # kJ/kg
    ho = hi - x*dhrefa - (1-x)*dhrefw
    return ho


def calcul_q_Tpx(Tg,p,x):
    REFPROP.set_mass_fractions([1 - x,x])
    def fct_Tk(q):
        REFPROP.update(PQ_INPUTS, p, q)
        Tk = REFPROP.T()
        return Tk - Tg
    qmole = brentq(fct_Tk,0,1)
    xdew = calcul_x_Tpd(Tg,p)
    xbubble = calcul_x_Tpb(Tg,p)
    Mv = xdew*Mammonia + (1-xdew)*Mwater
    Ml = xbubble*Mammonia + (1-xbubble)*Mwater
    rap = Ml/Mv
    qmass = qmole/(qmole+(1-qmole)*rap)
    return qmass



def calcul_x_Tpd(Tg,p):
    def fct(x):
        REFPROP.set_mass_fractions([1 - x,x])
        REFPROP.update(PQ_INPUTS, p, 1)
        Tk = REFPROP.T()
        return Tk - Tg
    z = brentq(fct,0,1)
    return z

def calcul_x_Tpb(Tg,p):
    def fct(x):
        REFPROP.set_mass_fractions([1 - x,x])
        REFPROP.update(PQ_INPUTS, p, 0)
        Tk = REFPROP.T()
        return Tk - Tg
    z = brentq(fct,0,1)
    return z
def calcul_x_Tpq(Tg,p,qmass):

    xdew = calcul_x_Tpd(Tg,p)
    xbubble = calcul_x_Tpb(Tg,p)
    Mv = xdew*Mammonia + (1-xdew)*Mwater
    Ml = xbubble*Mammonia + (1-xbubble)*Mwater
    rap = Mv/Ml
    qmole = qmass/(qmass+(1-qmass)*rap)
    def fct(x):
        REFPROP.set_mass_fractions([1 - x,x])
        REFPROP.update(PQ_INPUTS, p, qmole)
        Tk = REFPROP.T()
        return Tk - Tg
        xdew = calcul_x_Tpd(Tg,p)
    z = brentq(fct,0,1)
    return z

def PropsSI_nh3_h20(prop_out, prop_in1, arg1,prop_in2, arg2, x = 0.5):
    REFPROP.set_mass_fractions([1 - x,x])
    prop_in1.upper()
    prop_out.upper()
    prop_in2.upper()
    if prop_in1 == 'T'or prop_in2 == 'T':
        if prop_in2 == 'T':
            a = arg1
            arg1 = arg2
            arg2 = a
            prop_in2 = prop_in1
            prop_in1 = 'T'
        if prop_in2 == 'P':
            REFPROP.update(PT_INPUTS, arg2, arg1)
            if prop_out == 'H':
                y = REFPROP.hmass()
            if prop_out == 'S':
                y = REFPROP.smass()
            if prop_out == 'Q':
                y = REFPROP.Q()
        if prop_in2 == 'Q':
            REFPROP.update(QT_INPUTS, arg2, arg1)
            if prop_out == 'P':
                y = REFPROP.p()
            if prop_out == 'H':
                y = REFPROP.hmass()
            if prop_out == 'S':
                y = REFPROP.smass()
        if prop_in2 == 'H':
            REFPROP.update(HmassT_INPUTS, arg2, arg1)
            if prop_out == 'P':
                y = REFPROP.p()
            if prop_out == 'S':
                y = REFPROP.smass()
            if prop_out == 'Q':
                y = REFPROP.Q()
        if prop_in2 == 'S':
            REFPROP.update(SmassT_INPUTS, arg2, arg1)
            if prop_out == 'P':
                y = REFPROP.p()
            if prop_out == 'H':
                y = REFPROP.hmass()
            if prop_out == 'Q':
                y = REFPROP.Q()
    elif prop_in1 == 'P'or prop_in2 == 'P':
        if prop_in2 == 'P':
            a = arg1
            arg1 = arg2
            arg2 = a
            prop_in2 = prop_in1
            prop_in1 = 'P'
        if prop_in2 == 'Q':
            REFPROP.update(PQ_INPUTS, arg1, arg2)
            if prop_out == 'T':
                y = REFPROP.T()
            if prop_out == 'H':
                y = REFPROP.hmass()
            if prop_out == 'S':
                y = REFPROP.smass()
        if prop_in2 == 'H':
            REFPROP.update(HmassP_INPUTS, arg2, arg1)
            if prop_out == 'T':
                y = REFPROP.T()
            if prop_out == 'S':
                y = REFPROP.smass()
            if prop_out == 'Q':
                y = REFPROP.Q()
        if prop_in2 == 'S':
            REFPROP.update(PSmass_INPUTS, arg1, arg2)
            if prop_out == 'T':
                y = REFPROP.T()
            if prop_out == 'H':
                y = REFPROP.hmass()
            if prop_out == 'Q':
                y = REFPROP.Q()
    return y



