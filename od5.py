from odlib import *
import numpy as np
from math import sin, cos, sqrt, asin, acos

'''
Calculate the orbital elements of any object rotating around the Sun, given its right ascension and 
declination on three separate days (time in Julian days) and the equatorial Earth-Sun vectors on those days (in AU)
'''

# MATH CONSTANTS
k = 0.01720209847
c = 173.145


# AIJ DATA
RA = [RAdecimalToHMS(17.633096*15), RAdecimalToHMS(17.642617*15),
      RAdecimalToHMS(17.661737*15)]
DEC = [DECdecimalToDMS(-24.906182), DECdecimalToDMS(-22.779707),
       DECdecimalToDMS(-21.790782)]
t1 = 2458671.778
t2 = 2458679.772
t3 = 2458683.774
R = [[-2.568180989317604E-01,  9.026400560686667E-01,  3.912553385188193E-01],
       [-3.849247965056879E-01, 8.631971645323986E-01, 3.741557451122443E-01],
      [-4.465913217555757E-01, 8.375970194890799E-01,  3.630559556728792E-01]] #in equatorial coordinates REAL
expected_elements = [2.640226301080940 , 4.059008132819690E-01, 9.527069429767748,
                     2.809685769461632E+02 , 3.567195418376240E+02, 1.777935844299149] # EXPECTED VALUES OF a, e, i, big om, w, M


# Convert RA and DEC coordinates to radians
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


# Set starting values for approximation
r2 = [3]
r2dot = [3]
r2o = [0]
r2doto = [0]


# f series function
def f(tau):
    series_pt1 = 1-(tau**2)/(2*mag(r2)**3) + \
                 (tau**3)*dot(r2, r2dot)/(2*mag(r2)**5)
    coeff = (tau**4)/(24*mag(r2)**3)
    disgusting = 3*(dot(r2dot, r2dot)/(mag(r2)**2) - 1/(mag(r2)**3)) - \
                 15*(dot(r2, r2dot)/(mag(r2)**2))**2 + 1/(mag(r2)**3)
    return series_pt1 + coeff*disgusting

# g series function
def g(tau):
    return tau - (tau**3)/(6*(mag(r2)**3)) + (tau**4)*(dot(r2, r2dot))/(4*mag(r2)**5)


# Determine if two vectors are similar enough to be considered equivalent
def determineClose(v1, v2):
    for i in range(len(v1)):
        if(abs(v1[i]-v2[i]) > 10**-11):
            return False
    return True


# Method of Gauss algorithm to approximate r2 and r2dot
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
    count += 1

      
# Convert r2 vectors to the ecliptic plane
r2 = rotate_x(r2, 23.4367505323)
r2dot = rotate_x(r2dot, 23.4367505323)


# Calculate orbital elements based on r2
def calc_orbital_elements(r, rdot):
    vsquared = dot(rdot, rdot)
    h = cross(r, rdot)
    print(h)
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


# Format output of orbital elements
def format(expected, actual, errors):
    print("Semimajor axis (a):")
    print("Expected value:", expected[0], "Calculated value:", actual[0], \
          "Percent error:", errors[0]*100,"%")
    print()
    print("Eccentricity (e):")
    print("Expected value:", expected[1], "Calculated value:", actual[1], \
          "Percent error:", errors[1]*100,"%")
    print()
    print("Inclination (i):")
    print("Expected value:", expected[2], "Calculated value:", actual[2], \
          "Percent error:", errors[2]*100,"%")
    print()
    print("Longitude of ascending node (OM):")
    print("Expected value:", expected[3], "Calculated value:", actual[3], \
          "Percent error:", errors[3]*100,"%")
    print()
    print("Argument of perihelion (w):")
    print("Expected value:", expected[4], "Calculated value:", actual[4], \
          "Percent error:", errors[4]*100,"%")
    print()
    print("Mean anomaly (M):")
    print("Expected value:", expected[5], "Calculated value:", actual[5], \
          "Percent error:", errors[5]*100,"%")


def main(r, rdot, expected):
    elements = calc_orbital_elements(r, rdot)
    errors = calc_percent_error(expected, elements)
    format(expected, elements, errors)


main(r2, r2dot, expected_elements)
