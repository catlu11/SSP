from odlib import *
import numpy as np
from math import sin, cos, sqrt, asin, acos

# Math constants
k = 0.01720209847
c = 173.145

def f(tau, r2, r2dot):
    series_pt1 = 1-(tau**2)/(2*mag(r2)**3) + \
                 (tau**3)*dot(r2, r2dot)/(2*mag(r2)**5)
    coeff = (tau**4)/(24*mag(r2)**3)
    disgusting = 3*(dot(r2dot, r2dot)/(mag(r2)**2) - 1/(mag(r2)**3)) - \
                 15*(dot(r2, r2dot)/(mag(r2)**2))**2 + 1/(mag(r2)**3)
    return series_pt1 + coeff*disgusting
def g(tau, r2, r2dot):
    return tau - (tau**3)/(6*(mag(r2)**3)) + (tau**4)*(dot(r2, r2dot))/(4*mag(r2)**5)

def determineClose(v1, v2):
    for i in range(len(v1)):
        if(abs(v1[i]-v2[i]) > 10**-11):
            return False
    return True

def calc_r():
    RA = [RAdecimalToHMS(np.random.normal(17.633096*15, 0.2675/3600)),
          RAdecimalToHMS(np.random.normal(17.642617*15, 0.2675/3600)),
          RAdecimalToHMS(np.random.normal(17.661737*15, 0.2675/3600))]
    DEC = [DECdecimalToDMS(np.random.normal(-24.906182, 0.2960/3600)),
           DECdecimalToDMS(np.random.normal(-22.779707, 0.2960/3600)),
           DECdecimalToDMS(np.random.normal(-21.790782, 0.2960/3600))]
    t1 = 2458671.778
    t2 = 2458679.772
    t3 = 2458683.774
    R = [[-2.568180989317604E-01,  9.026400560686667E-01,  3.912553385188193E-01],
       [-3.849247965056879E-01, 8.631971645323986E-01, 3.741557451122443E-01],
      [-4.465913217555757E-01, 8.375970194890799E-01,  3.630559556728792E-01]]

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

    r2 = [3]
    r2dot = [3]
    r2o = [0]
    r2doto = [0]

    while(determineClose(r2, r2o) == False and determineClose(r2dot, r2doto) == False and count < 10000):
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
            r2dot = (rvec[2]-rvec[0])/tau0

        # Calculate f and g series
        f1 = f(tau1, r2, r2dot)
        f3 = f(tau3, r2, r2dot)
        g1 = g(tau1, r2, r2dot)
        g3 = g(tau3, r2, r2dot)

        # Update a and b
        a1 = g3/(f1*g3-f3*g1)
        a3 = g1/(f3*g1-f1*g3)
        b1 = f3/(f3*g1-f1*g3)
        b3 = f1/(f1*g3-f3*g1)

        # Update r2dot
        r2doto = r2dot
        r2dot = b1*rvec[0] + b3*rvec[2]
        count += 1

    r2 = rotate_x(r2, 23.4367505323)
    r2dot = rotate_x(r2dot, 23.4367505323)
    return [r2, r2dot]

def calc_orbital_elements(r, rdot):
    vsquared = dot(rdot, rdot)
    h = cross(r, rdot)

    # Calculate semimajor axis
    a = 1/(2/mag(r) - vsquared)

    # Calculate eccentricity
    ec = sqrt(1-(mag(h)**2)/a)

    # Calculate inclination
    i = acos(h[2]/mag(h))
    
    # Calculate longitude of the ascending node
    sin_omega = h[0]/(mag(h)*sin(i))
    cos_omega = -h[1]/(mag(h)*sin(i))
    l_omega = trig_to_rad(sin_omega, cos_omega)
    if(l_omega < 0):
        l_omega += (2*pi)    
        
    # Calculate argument of perihelion
    sin_nu = (a*(1-ec**2)/mag(h))*(dot(r, rdot)/mag(r))*(1/ec)
    cos_nu = ((a*(1-ec**2)/mag(r))-1)*(1/ec)
    sin_u = r[2]/(mag(r)*sin(i))
    cos_u = (r[0]*cos_omega + r[1]*sin_omega)/mag(r)
    u = trig_to_rad(sin_u, cos_u)
    nu = trig_to_rad(sin_nu, cos_nu)
    p_omega = u - nu
    if(p_omega < 0):
        p_omega += (2*pi)
        
    # Calculate mean anomaly
    E = acos((1/ec)*(1-mag(r)/a))
    M = E - ec*sin(E)

    return [a, ec, i*180/pi, l_omega*180/pi, p_omega*180/pi, M*180/pi]

a = []
ec = []
inc = []
om = []
w = []
M = []
def main():
    for i in range(20):
        rvs = calc_r()
        for i in range(3):
            rvs[0][i] = str(rvs[0][i])
            rvs[1][i] = str(rvs[1][i]*k*365.25)    
        print(' '.join(list(rvs[0])))
        print(' '.join(list(rvs[1])))
        print()
##        o = calc_orbital_elements(rvs[0], rvs[1])
##        a.append(o[0])
##        ec.append(o[1])
##        inc.append(o[2])
##        om.append(o[3])
##        w.append(o[4])
##        M.append(o[5])

main()
##elements = [a, ec, inc, om, w, M]
##print(mean(a), mean(ec), mean(inc), mean(om), mean(w), mean(M))
##np.savetxt('elements.csv', elements, delimiter = ',')
