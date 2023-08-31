# coding utf-8
from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse
from hydraulic_md import *


f = open("pac6m.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p6 = array([float(x) for x in ligne2])/9.8
deb6 = L1[0]
deb26 = 1.25*deb6
hpi6 = L1[1]
Wp6 = L1[2]
print('Wp 6= ',Wp6)

f = open("pac5m.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p5 = array([float(x) for x in ligne2])/9.8
deb5 = L1[0]
deb25 = 1.25*deb5
hpi5 = L1[1]
Wp5 = L1[2]
print('Wp 5= ',Wp5)
f = open("pac5bm.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p5b = [float(x) for x in ligne2]
deb5b = L1[0]
deb25b = 1.25*deb5b
hpi5b = L1[1]
Wp5b = L1[2]

print('Wp 5b = ',Wp5b)
#
#
f = open("pac4m.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p4 = array([float(x) for x in ligne2])/9.8
deb4 = L1[0]
deb24 = 1.25*deb4
hpi4 = L1[1]
Wp4 = L1[2]
print('Wp 4= ',Wp4)
#
f = open("pac4bm.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p4b = array([float(x) for x in ligne2])/9.8
deb4b = L1[0]
deb24b = 1.25*deb4b
hpi4b = L1[1]
Wp4b = L1[2]
print('Wp = 4b',Wp4b)
#
#
f = open("pac3m.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p3 = array([float(x) for x in ligne2])/9.8
deb3 = L1[0]
deb34 = 1.25*deb3
hpi3 = L1[1]
Wp3  = L1[2]
print('Wp = 3',Wp3)
#
f = open("pac3bm.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p3b = array([float(x) for x in ligne2])/9.8
deb3b = L1[0]
deb34b = 1.25*deb3b
hpi3b = L1[1]
Wp3b = L1[2]

print('Wp = 3b',Wp3b)
#
f = open("pac2m.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p2 = array([float(x) for x in ligne2])/9.8
deb2 = L1[0]
deb22 = 1.25*deb2
hpi2 = L1[1]
Wp2 = L1[2]
print('Wp = 2',Wp2)
#
f = open("pac2bm.txt", "r")
M = f.readlines()
ligne1 = M[0].rstrip().split(',')
L1 = [float(x) for x in ligne1]
ligne2 = M[1].rstrip().split(',')
p2b = array([float(x) for x in ligne2])/9.8
deb2b = L1[0]
deb22b = 1.25*deb2b
hpi2b = L1[1]
Wp2b = L1[2]
print('Wp = 2b',Wp2b)
x1 = array([2,3,4,5,6])
y1 = array([Wp2,Wp3,Wp4,Wp5,Wp6])
y2 = array([Wp2b,Wp3b,Wp4b,Wp5b,Wp6])
#
rcParams.update({'figure.autolayout': True})
y3 = ones(5)
y5 = zeros(5)
X = 0.52723
Z =  fct_lemire(1,X)
rend = array([0.48,0.55,0.58,0.59,0.6])
for i in range(0,5):
    q = (2 + i)/6
    y3[i] = fct_lemire(q,X)
    DH = 12.2*((1-X)*q**2 + X)
    y5[i] = q*4.8*9.81*DH/rend[i]
y4 = y3*y1[4]/Z
p1 = plot(x1,y4,'k-+',x1,y2,'k-o',x1,y5,'kx',markersize=13)
ax = gca()
dx = ax.get_xlim()
dy = ax.get_ylim()
ratio = (dy[1]-dy[0])/(dx[1] - dx[0])
rr = 0.05
grid(True,which='both')
#title('Coubes de pertes de charge  ')
fts = 16
ftst = 14
legend(('Lemire correlation','Exemple 8.1','Sfeir relation'),fontsize = fts)
xlabel('Number of Heat pumps',fontsize = fts)
ylabel(r'Power (W)',fontsize = fts)
xticks(fontsize=ftst)
yticks(fontsize=ftst)
ax = gca()
ax.set_xticks([2,3,4,5,6])
ax.set_xticklabels(['two','three','four','five','six'])
show()
