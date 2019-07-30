from odlib import *
from math import *
# OD code 1 - Orbital elements

def calc_orbital_elements(r, rdot):
    vsquared = dot(rdot, rdot)
    h = cross(r, rdot)

    #Calculate semimajor axis
    a = 1/((2/mag(r)) - vsquared)

    #Calculate eccentricity
    ec = sqrt(1-((mag(h)**2)/a))

    #Calculate inclination
    i = acos(h[2]/mag(h)) 

    #Calculate longitude of the ascending node
    sin_omega = h[0]/(mag(h)*sin(i))
    cos_omega = -h[1]/(mag(h)*sin(i))
    l_omega = trig_to_rad(sin_omega, cos_omega)
    if(l_omega < 0):
        l_omega += 2*pi
        
    #Calculate argument of perihelion
    sin_nu = (a*(1-ec**2)/mag(h))*(dot(r, rdot)/mag(r))*(1/ec)
    cos_nu = ((a*(1-ec**2)/mag(r))-1)*(1/ec)
    sin_u = r[2]/(mag(r)*sin(i))
    cos_u = (r[0]*cos_omega + r[1]*sin_omega)/mag(r)
    u = trig_to_rad(sin_u, cos_u)
    nu = trig_to_rad(sin_nu, cos_nu)
    p_omega = u - nu
    if(p_omega < 0):
        p_omega += 2*pi
        
    #Calculate mean anomaly
    E = acos((1/ec)*(1-mag(r)/a))
    M = E - ec*sin(E)

    return [a, ec, i*180/pi, l_omega*180/pi, p_omega*180/pi, M]

def calc_percent_errors(expected, actual):
    errors = []
    for i in range(len(expected)):
        errors.append(abs((actual[i]-expected[i])/expected[i]))
    return errors

def format(outputs, expected, errors):
    print("Semimajor axis (a):")
    print("Expected value:", expected[0], "Calculated value:", outputs[0], \
          "Percent error:", errors[0]*100)
    print()
    print("Eccentricity (e):")
    print("Expected value:", expected[1], "Calculated value:", outputs[1], \
          "Percent error:", errors[1]*100)
    print()
    print("Inclination angle (i):")
    print("Expected value:", expected[2], "Calculated value:", outputs[2], \
          "Percent error:", errors[2]*100)
    print()
    print("Longitude of ascending node (Q):")
    print("Expected value:", expected[3], "Calculated value:", outputs[3], \
          "Percent error:", errors[3]*100)
    print()
    print("Argument of perihelion (w):")
    print("Expected value:", expected[4], "Calculated value:", outputs[4], \
          "Percent error:", errors[4]*100)
    print()
    print("Mean anomaly (M):")
    print("Expected value:", expected[5], "Calculated value:", outputs[5], \
          "Percent error:", errors[5]*100)

def main(r, rdot, expected):
    elements = calc_orbital_elements(r, rdot)
    errors = calc_percent_errors(expected, elements)
    format(elements, expected, errors)
    
r = [0.3975379996283998, -1.224961944120962, 0.4748988248594287]
rdot = [0.0004601682762, 0.0001081844755, 0.0001514210609]
expected = [1.056, 0.344233, 25.155, 236.23, 255.5, 140.4]
main(r, rdot, expected)
