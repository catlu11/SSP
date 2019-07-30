from odlib import HMStoDeg, DMStoDeg, rotate_y, rotate_z, mag, trig_to_rad
import numpy as np
from math import sin, cos, asin, pi
# LST, RA, DEC 
RA = [[20,0,0], [19,42,0], [18,10,0], [19,57,0], [18,42,0], [19,5,0], [22,58,0]]
DEC = [[-26,0,0], [-21,0,0], [-22,30,0], [2,0,0], [-17,0,0], [-20,0,0],[-18,0,0]]
LST = [286.5, 277.5, 232.5, 292.5, 273.75, 241, 297]

# OD3: RA and DEC to x, y, z
RA_rad = [HMStoDeg(x[0], x[1], x[2]) * pi/180 for x in RA]
DEC_rad = [DMStoDeg(x[0], x[1], x[2]) * pi/180 for x in DEC]

rho_x = [cos(RA_rad[x])*cos(DEC_rad[x]) for x in range(len(RA_rad))]
rho_y = [sin(RA_rad[x])*cos(DEC_rad[x]) for x in range(len(RA_rad))]
rho_z = [sin(DEC_rad[x]) for x in range(len(RA_rad))]
rho = [[rho_x[x], rho_y[x], rho_z[x]] for x in range(len(RA))]

# Rotate about z axis CCW by LST then about y axis CW by 90-40=50 degrees
rho_rot1 = []
rho_rot2 = []
for i in range(len(LST)):
    rho[i] = rotate_y(rotate_z(rho[i], -LST[i]), 50)

# Calculate alt
alt = [asin(x[2])*180/pi for x in rho]

# Calculate azimuth
azi = []
for i in range(len(alt)):
    azi.append(trig_to_rad(rho[i][1]/cos(alt[i]*pi/180), rho[i][0]/cos(alt[i]*pi/180))*180/pi)
print(azi)
