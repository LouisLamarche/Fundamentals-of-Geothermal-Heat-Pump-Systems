#
#
from numpy import *
from geothermal_mod import *
from conversion_mod import *
from  matplotlib.pyplot import *
from  scipy.optimize import newton
#
pcl= array([0.006,-0.42,11.5])
q_cl = 20000.
Tci = 14.0   #
rho = 1000
mph = 0.8   # debit  kg/s si o a lde l'eau)
Ch = mph*4200.          # en kW/K
# climatisation
UA = 4000.0
mpcv = array([0.6,0.7,0.8,0.9])
nn = len(mpcv)
Wtot = zeros(nn)
for i in range(0,nn):
    mpc = mpcv[i]
    Cc = mpc*4200
    Cmin = min(Cc,Ch)
    Cmax = max(Cc,Ch)
    Cr = Cmin/Cmax
    NTU = UA/Cmin
    if Cr < 1:
        y = exp(-NTU*(1-Cr))
        ep = (1-y)/(1-Cr*y)
    else:
        ep = NTU/(1+NTU)

    def Calcul_Thi(Thi):
        qmax = Cmin*(Thi-Tci)
        q = ep*qmax
        Tho = Thi - q/Ch
        coph = polyval(pcl,Tho)
        qech = q_cl*(coph+1)/coph
        y = qech - q
        return y

    Thi = newton(Calcul_Thi,18)
    qmax = Cmin*(Thi-Tci)
    q = ep*qmax
    Tho = Thi - q/Ch
    Tco = Tci + q/Cc
    copcl = polyval(pcl,Tho)
    qech = q_cl*(copcl+1)/copcl
    Wcomp =q_cl/copcl
    beta = 1000.0   # venant de 11.3
    C= 2 # min2/m5
    C = 2*3600  # sec2/m5
    Ccond = 15e6
    Q = mpc/rho  # m3/s
    st = beta*Q + C*Q**2
    ho = 10.0
    hp = Ccond*Q**2
    ht = ho + st + hp
    Wf = mpc*9.8*ht
    rend = 0.5
    Wpompe = Wf/rend
    print ('%0.0f' %  (Wcomp+Wpompe))
    print ('%0.1f' %  Tco)
    Wtot[i] = Wcomp+Wpompe
    COPsys = q_cl/(Wcomp+Wpompe);
plot(mpcv,Wtot)
show()