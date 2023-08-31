#coding: utf-8
#
import numpy as np

def pwf(N,i,d,cas = 'end'):
    if (cas == 'end'):
        if i == d :
            vaf = N/(1.0+d)
        else:
            vaf = 1.0/(d-i)*(1.0-((1.0+i)/(1.0+d))**N)
    else:
        if i == d :
            vaf = N
        else:
            vaf = (1.0+i)/(d-i)*(1.0-((1.0+i)/(1.0+d))**N)
    return vaf



def fwf(N,i,d,cas = 'end'):
    vaf = pwf(N,i,d,cas)*(1+d)**N
    return vaf


def pwf_int(Np,Na,m,d,cas = 'end'):
    Nmin = min(Np,Na)
    x1 = pwf(Nmin,0,d,cas)
    x2 = pwf(Np,0,m,cas)
    x3 = pwf(Nmin,m,d,cas)
    vaf = (x1/x2+x3*(m-1/x2))
    return vaf

def present_value(A,N,ii,t,cas = 'end'):
    if(cas == 'end'):
        va = A*(1.0+ii)**(N-1)/(1+t)**N
    else:
        va = A*(1.0+ii)**(N)/(1+t)**N
    return va
def future_value(A,N,t):
    va = A*(1+t)**N
    return va




def mirr(dpa,gain,rr,fr):
    num= 0
    den = dpa
    n = len(gain)
    for i in range(0,n):
        if gain[i] > 0:
            num = num + gain[i]*(1 + rr)**(n-1-i)
        else:
            den = den - gain[i]/(1 + fr)**(i+1)
    mirr = (num/den)**(1/n) -1
    return mirr