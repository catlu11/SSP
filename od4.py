from odlib import HMStoDeg, DMStoDeg, mag, dot, cross, rotate_x, mean
import numpy as np
from math import *
# OD 4

# Math constants
k = 0.01720209847
c = 173.145

# JPL data
RA = [[18,41,48.08],[18,14,30.53], [17,54,29.36]]
DEC = [[36,15,14.0],[36,9,46], [34,29,32.2]]
R = [[-0.2069177625542078, 0.9132792377095497, 0.3959074080096284],
      [-0.3689266994628273, 0.8690904003976556, 0.3767540914983889],
     [-0.5205074019556455, 0.8003914633545058, 0.3469760777959491]]
t1 = 2458303.5
t2 = 2458313.5
t3 = 2458323.5

# Convert data to radians
RA_rad = [HMStoDeg(x[0], x[1], x[2]) * pi/180 for x in RA]
DEC_rad = [DMStoDeg(x[0], x[1], x[2]) * pi/180 for x in DEC]

# Calculate rho hat from RA and DEC values
rho_x = [cos(RA_rad[x])*cos(DEC_rad[x]) for x in range(len(RA_rad))]
rho_y = [sin(RA_rad[x])*cos(DEC_rad[x]) for x in range(len(RA_rad))]
rho_z = [sin(DEC_rad[x]) for x in range(len(RA_rad))]
rho_hat_calc = [[rho_x[x], rho_y[x], rho_z[x]] for x in range(len(RA))]
rho_hat = np.array(rho_hat_calc)

# Adjust for time difference w/ taus
tau1 = k*(t1-t2)
tau3 = k*(t3-t2)
tau0 = tau3 - tau1

# Approximate a1 and a3
a1 = tau3/tau0
a3 = -tau1/tau0
count = 0

r2 = [5]
r2dot = [5]
r2o = [0]
r2doto = [0]

def f(tau):
    series_pt1 = 1-(tau**2)/(2*mag(r2)**3) + \
                 (tau**3)*dot(r2, r2dot)/(2*mag(r2)**5)
    coeff = (tau**4)/(24*mag(r2)**3)
    disgusting = 3*(dot(r2dot, r2dot)/(mag(r2)**2) - 1/(mag(r2)**3)) - \
                 15*(dot(r2, r2dot)/(mag(r2)**2))**2 + 1/(mag(r2)**3)
    print(coeff*disgusting)
    return series_pt1 + coeff*disgusting
def g(tau):
    return tau - (tau**3)/(6*(mag(r2)**3)) + (tau**4)*(dot(r2, r2dot))/(4*mag(r2)**5)

def determineClose(v1, v2):
    for i in range(len(v1)):
        if(abs(v1[i]-v2[i]) > 10**-11):
            return False
    return True

while(determineClose(r2, r2o) == False and determineClose(r2dot, r2doto) == False):
    D0 = dot(rho_hat[0], cross(rho_hat[1], rho_hat[2]))
    D = np.ones((3,3))
    for j in range(3):
        D[0,j] = dot(cross(R[j], rho_hat[1]), rho_hat[2])
        D[1,j] = dot(cross(rho_hat[0], R[j]), rho_hat[2])
        D[2,j] = dot(rho_hat[0], cross(rho_hat[1], R[j]))
    
    # Calculate rhos
    rho = np.ones(3)
    rho[0] = (a1*D[0,0]+(-1)*D[0,1]+a3*D[0,2])/(a1*D0)
    rho[1] = (a1*D[1,0]+(-1)*D[1,1]+a3*D[1,2])/(-1*D0)
    rho[2] = (a1*D[2,0]+(-1)*D[2,1]+a3*D[2,2])/(a3*D0)

    # Adjust taus
    tau1 = k*((t1-rho[0]/c)-(t2-rho[1]/c))
    tau3 = k*((t3-rho[2]/c)-(t2-rho[1]/c))
    tau0 = tau3 - tau1
    
    # Approximate vectors r1, r2, r3, r2dot
    rvec = []
    for i in range(3):
        rvec.append(rho[i]*rho_hat[i] - R[i])
    rvec = np.array(rvec)
    r2o = r2
    r2 = rvec[1]
    if(count == 0):
        count += 1
        r2dot = (rvec[2]-rvec[0])/tau0

    # Calculate f and g series
    f1 = f(tau1)
    f3 = f(tau3)
    g1 = g(tau1)
    g3 = g(tau3)

    # Update a and b
    a1 = g3/(f1*g3-f3*g1)
    a3 = g1/(f3*g1-f1*g3)
    b1 = f3/(f3*g1-f1*g3)
    b3 = f1/(f1*g3-f3*g1)

    # Update r2dot
    r2doto = r2dot
    r2dot = b1*rvec[0] + b3*rvec[2]


r2 = rotate_x(r2, 23.4367505323)
r2dot = rotate_x(r2dot, 23.4367505323)
er2 = [0.3970630876567, -1.22507372544802, 0.4747425864661720]
er2dot = [0.662642317731134, 0.155785644680128, 0.218046327693182]

xdif = abs(r2[0]-er2[0]) / mean([r2[0], er2[0]])
ydif = abs(r2[1]-er2[1]) / mean([r2[1], er2[1]])
zdif = abs(r2[2]-er2[2]) / mean([r2[2], er2[2]])
errors = np.array([xdif, ydif, zdif])

xddif = abs(r2dot[0]-er2dot[0]) / mean([r2dot[0], er2dot[0]])
yddif = abs(r2dot[1]-er2dot[1]) / mean([r2dot[1], er2dot[1]])
zddif = abs(r2dot[2]-er2dot[2]) / mean([r2dot[2], er2dot[2]])
errorsd = np.array([xddif, yddif, zddif])

print("Expected r:", er2)
print()
print("Calculated r:", r2)
print()
print("Percent difference:", errors*100)
print()
print("Expected rdot:", er2dot)
print()
print("Calculated rdot:", r2dot)
print()
print("Percent difference:", errorsd*100)
print()

