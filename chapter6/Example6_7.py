#coding: utf-8
#
# Exemple 6.7 written by Louis lamarche 29 septemeber 2017
#
from geothermal_md import *
import numpy as np
#
# data
#
als = 0.1 # m2/jr
ks = 2.0
alhr = als/24.0
qp = 8.0
d = 6.0
n_years = 10
H = 100.0
nx = 3
ny = 3
Tp1 = Tp_ashrae(nx,ny,d,qp,als,ks,n_years)
Tp1b = Tp_Bernier(nx,ny,d,qp,als,ks,n_years,H)        # A = 1
Tp1c = Tp_Capozza(nx,ny,d,qp,als,ks,n_years,H=H)
Tp1f = Tp_eight(nx,ny,d,qp,als,ks,n_years,H=H,conf = 'R')
nx = 9
ny = 1
Tp2 = Tp_ashrae(nx,ny,d,qp,als,ks,n_years)
Tp2b = Tp_Bernier(nx,ny,d,qp,als,ks,n_years,H)        # A = 1
Tp2c = Tp_Capozza(nx,ny,d,qp,als,ks,n_years,H=H)
Tp2f = Tp_eight(nx,ny,d,qp,als,ks,n_years,H=H,conf = 'NR')
print ('Tp 3x3 (Philippe) = {:.2f}'.format(Tp1b))
print ('Tp 3x3 (Capoza) = {:.2f}'.format(Tp1c))
print ('Tp 3x3 (Fossa) = {:.2f}'.format(Tp1f))
print ('Tp 9x1 (Philippe) = {:.2f}'.format(Tp2b))
print ('Tp 9x1 (Capozza) = {:.2f}'.format(Tp2c))
print ('Tp 9x1 (Fossa) = {:.2f}'.format(Tp2f))