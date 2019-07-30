from math import sqrt, sin, cos
from odlib import *
import numpy as np

# a, e, i, big omega, lil omega, M
orbital_elements = [2.640226301080940 , 4.059008132819690E-01, 9.527069429767748,
                     2.809685769461632E+02 , 3.567195418376240E+02, 1.777935844299149]
R = [-4.465913217555757E-01, 8.375970194890799E-01,  3.630559556728792E-01]
expected = [[17,39,42.29], [-21,47,27.5]]
t = 2458683.774

def f(x, e, M):
    return x-e*sin(x)-M
def fprime(x, e):
    return 1-e*cos(x)

def generate_coords(orbital_elements, t):
    k = 0.01720209847
    t2 = 2458679.772
    n = sqrt(1/(orbital_elements[0])**3)
    T = t2 - (orbital_elements[5]*pi/180)/(n*k)

    # Calculate M
    M = n*k*(t-T)
    M = M % (2*pi)
    
    # Calculate E
    e = orbital_elements[1]
    E = 0
    while(abs(f(E, e, M)) > (10**(-15))):
        fx = f(E, e, M)
        d = fprime(E, e)
        E = E-fx/d
    # Calculate r
    r = [orbital_elements[0]*cos(E) - orbital_elements[0]*orbital_elements[1],
         orbital_elements[0]*sqrt(1-orbital_elements[1]**2)*sin(E),
         0]

    # First rotation
    r = rotate_z(r, orbital_elements[4])

    # Second rotation
    r = rotate_x(r, -orbital_elements[2])

    # Third rotation
    r = rotate_z(r, orbital_elements[3])

    ## Fourth rotation
##    r = rotate_x(r, -23.4367505323)

    # Calculate rho and rho hat
##    rho = np.array(R) + np.array(r)
##    rho_hat = rho/mag(rho)

    # Calculate RA and DEC
##    DEC = asin(rho_hat[2])
##    RA = trig_to_rad(rho_hat[1]/cos(DEC), rho_hat[0]/cos(DEC))

    return r

##RA = (RA*180/pi)
##DEC = (DEC*180/pi)

# Calculate percent errors
##RAdif = abs(HMStoDeg(expected[0][0], expected[0][1], expected[0][2])-RA) \
##        /mean([HMStoDeg(expected[0][0], expected[0][1], expected[0][2]),RA])
##DECdif = abs(DMStoDeg(expected[1][0], expected[1][1], expected[1][2])-DEC) \
##        /mean([DMStoDeg(expected[1][0], expected[1][1], expected[1][2]),DEC])
##
##print("Calculated: RA =", RA," DEC =", DEC)
##print("Expected: RA =", HMStoDeg(expected[0][0], expected[0][1], expected[0][2])," DEC =",
##DMStoDeg(expected[1][0], expected[1][1], expected[1][2]))
##print("RA error = "+str(RAdif*100)+"%  DEC error = "+str(DECdif*100)+"%")
