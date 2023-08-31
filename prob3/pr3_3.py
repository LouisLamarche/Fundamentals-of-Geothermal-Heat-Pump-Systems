Tinf = 0
UA12 = .4
UA1 = 0.6
qin1 = 0
qin2 = 12
T1 = 22
T2 = 22
COP = 3
#
# zone 2
#
qevap = qin2 + UA12*(T1 - T2)
#
# zone 1
qout = UA1*(T1 - Tinf)
qcond = COP*qevap/(COP-1)
qaux = qout - qcond
if qaux > 0:
    print('auxillary heat = ', qaux,' kW')
else:
    print('rejected  heat = ', -qaux,' kW')
#



