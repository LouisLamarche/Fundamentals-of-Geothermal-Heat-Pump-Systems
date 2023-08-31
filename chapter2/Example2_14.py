import numpy as np
# exemple 2.14
data = np.loadtxt("..\\data\\Temp_Mtl.txt")
Tamb1 = data[:,1]    # températures extérieures
hrm = np.array([744,672,744,720,744,720,744,744,720,744,720,744])
hrr = np.array([744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760])
dh = np.zeros(12)
i1 = 0
Ktot =  31.0
gains = 225.0
Ti = 22.0
Te = Ti - gains/Ktot
for im in range (0,12):
    for i in range (i1,hrr[im]):
        dh[im] = dh[im] + max(Te-Tamb1[i],0)
    i1 = hrr[im] + 1
q = dh*Ktot/1000.0   # vakeurs en kWh
dj = dh/24.0
dja = sum(dj)
qt = sum(q)
for im in range (0,12):
    ij = im + 1
    print ('Degree-days of month ' + '%d' % ij + ' are  :  ' + '%d' % dj[im])
for im in range (0,12):
    ij = im + 1
    print ('Heat losses  of month ' + '%d' % ij + ' are  :  ' + '%.2f' % q[im] + ' kWh')
print ('c) Total heat losses are  ')
print ('%.2f' % qt + ' kWh')