
from geothermal_md import *
from  matplotlib.pyplot import *

import numpy as np

cas = 'a'
if cas == 'a':
    Fo = 0.5
    tit = 'Fo = 0.5'
else:
    Fo = 50
    tit = 'Fo = 50'
r = np.arange(0.1,2.1,0.1)
n = len(r)
rr = np.arange(1,2.1,0.1)
n2 = len(rr)
y1 = np.zeros(n)
y2 = np.zeros(n2)
y3 = np.zeros(n2)
for i in range(0,n):
    y1[i] = G_function_isc(Fo,r[i])
for i in range(0,n2):
    y2[i] = G_function_ils(Fo,rr[i])
    y3[i] = G_function_ics(Fo,rr[i])

p1 = plot(r,y1,color = 'k',linestyle = '-',marker = 'o',markersize=8,label = 'ISC')
p2 = plot(rr,y2,color = 'k',linestyle = '-',marker = '+',markersize=8,label = 'ILS')
p3 = plot(rr,y3,color = 'k',linestyle = '-',marker = 'x',markersize=8,label = 'ICS')
legend(fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14

xlabel(' $\\tilde{r}$',fontname='Times new Roman',fontsize = fts)
ylabel(' G',fontsize = fts,fontname='Times new Roman')
title(tit,fontsize = 22,fontname='Times new Roman')

xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()
