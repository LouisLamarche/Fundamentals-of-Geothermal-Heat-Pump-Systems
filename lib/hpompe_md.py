import  numpy as np
def pompe_2007():   # q est en m3/s, h est en m
    p = np.array([-4.09482005e+07, -2.30170131e+05,  4.06869048e+00,  7.69450511e+01])
    return p
def h_pompe_1121(q):   # q est en l/s, h est en m
    p =  np.array([-0.06885822, -0.70165262, -3.24241969,  6.72386013])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def pompe_1121():   # q est en m3/s, gh est en m
    p =  np.array([-6.75499122e+08, -6.88321220e+06, -3.18081371e+04,  6.59610679e+01])
    return p
def gh_pompe_1121(q):   # q est en m3/s, gh est en m
    p =  np.array([-6.75499122e+08, -6.88321220e+06, -3.18081371e+04,  6.59610679e+01])
    h =  np.polyval(p,q)
    return h
def gh_pompe_0011(q):   # q est en m3/s, gh est en m
    p =  np.array([-6.59447711e+08, -3.82705350e+06, -3.82977535e+04,  9.47010593e+01])
    h =  np.polyval(p,q)
    return h
def pompe_0011():   # q est en m3/s, gh est en m
    p =  np.array([-6.59447711e+08, -3.82705350e+06, -3.82977535e+04,  9.47010593e+01])
    return p
def h_pompe_0014(q):   # q est en l/s, h est en m
    p =  np.array([-6.76191616e-03, -2.09580266e-01, -2.98742960e+00,  6.95780616e+00])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def pompe_0014():   # q est en m3/s, gh est en m
    p =  np.array([-6.63343975e+07, -2.05598241e+06, -2.93066844e+04,  6.82560784e+01])
    return p
def gh_pompe_0014(q):   # q est en m3/s, gh est en m
    p =  np.array([-6.63343975e+07, -2.05598241e+06, -2.93066844e+04,  6.82560784e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_1111(q):   # q est en l/s, h est en m
    p =  np.array([-0.00922663, -1.20814489, -0.89491534,  2.78492184])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def gh_pompe_1111(q):   # q est en m3/s, gh est en m
    p =  np.array([-9.05132874e+07, -1.18519014e+07, -8.77911946e+03,  2.73200833e+01])
    h =  np.polyval(p,q)
    return h
def h_pompe_1911_4_25(q):   # q est en l/s, h est en m
    p =  np.array([-0.01587967, -0.15556453, -0.04935972,  5.49522399])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_E90_2AAC_4_875(q):   # q est en l/s, h est en m
    p =  np.array([-2.04738735e-03,  -4.91155304e-02,   1.86622111e-01,    2.88403425e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_15101_5AD_6_25(q):   # q est en l/s, h est en m
    p =  np.array([-2.58713339e-03,  -3.40275632e-02,   9.20223405e-02,         1.25418762e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_1207(q):   # q est en l/s, h est en m
    p =  np.array([-7.95544363e-03, -1.77260106e-02, -9.64557475e-02,  1.39221865e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def gh_pompe_1207(q):   # q est en m3/s, gh est en m
    p =  np.array([-7.80429020e+07, -1.73892164e+05, -9.46230883e+02,  1.36576649e+02])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def pompe_1207():   # q est en m3/s, gh est en m
    p =  np.array([-7.80429020e+07, -1.73892164e+05, -9.46230883e+02,  1.36576649e+02])
    return p


def h_pompe_15102BD_7_875(q):   # q est en l/s, h est en m
    p =  np.array([-2.68499683e-03,  -1.18806380e-02,   7.52461582e-02,         2.00342024e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_15101_5BC_9_25(q):   # q est en l/s, h est en m
    p =  np.array([-1.19152324e-02,   4.98865190e-02,  -8.48601049e-02,3.10839719e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_15102g_12_5(q):   # q est en l/s, h est en m
    p =  np.array([-4.39109570e-03,   5.13923262e-02,  -3.09868224e-01, 5.15724980e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_15102g_12_5b(q):   # q est en l/s, h est en m
    p =  np.array([-3.62100539e-03,   2.73568106e-02,  -1.33474861e-01,
         5.14083217e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def rend_pompe_15102g_12_5(q):   # q est en l/s, h est en m
    p =  np.array([-4.39532773e-03,  -1.42390670e-01,   6.18653913e+00,
         1.53195015e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h

def h_pompe_1510_5_5(q):   # q est en l/s, h est en m
    p =  np.array([-3.87841882e-03,  -3.88589649e-02,1.31251196e-01,8.85130440e+00])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_15102bc_7_5(q):   # q est en l/s, h est en m
    p =   np.array([-3.00169007e-03,  -9.16273346e-03,   1.00398155e-01,1.73590368e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h

def h_pompe_1915_5_3(q):   # q est en l/s, h est en m
    p =  np.array([-0.22399274,0.34680112,-0.17705912,8.29803125])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_1915_5_75(q):   # q est en l/s, h est en m
    p =  np.array([-1.32401769e-01,1.14917400e-01,5.66071525e-03,1.01904275e+01])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def h_pompe_1915(q):   # q est en l/s, h est en m
    p =  np.array([-0.07371527, -0.0822366 , -0.05395484, 12.33191335])
    h =  np.polyval(p,q)                            # h est en metres
    return h
def gh_pompe_1915(q):   # q est en m3/s, gh est en m
    p =  np.array([-7.23146825e+08, -8.06741017e+05, -5.29296992e+02,  1.20976070e+02])
    h =  np.polyval(p,q)                            # h est en metres
    return h

def h_pompe_1935(q):   # q est en l/s, h est en m
    p = np.array([-0.07371527, -0.0822366 , -0.05395484, 12.33191335])
    h = np.polyval(p,q)                            # h est en metres
    return h
def gh_pompe_1935(q):   # q est en m3/s, gh est en m
    p = np.array([-1.21547184e+07, -9.03212269e+05,  1.83702735e+03,  6.83208441e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h
def pompe_1935():   # q est en m3/s, gh est en m
    p = np.array([-1.21547184e+07, -9.03212269e+05,  1.83702735e+03,  6.83208441e+01])
    return p
def gh_pompe_1935(q):   # q est en m3/s, gh est en m
    p = np.array([-1.21547184e+07, -9.03212269e+05,  1.83702735e+03,  6.83208441e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h
def pompe_1935b():   # q est en m3/s, gh est en m
    p = np.array([-1.36132846e+07, -1.01159774e+06,  2.05747063e+03,  7.65193454e+01])
    return p
def h_pompe_e90(q):   # q est en l/s, h est en m
    h = -5.74157895e-03*q**3 -1.79050043e-01*q**2 + 2.17708285e-01*q +  2.47621048e+01
    return h
