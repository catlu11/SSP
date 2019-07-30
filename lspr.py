from math import *
from copy import deepcopy
from odlib import *
import matplotlib.pyplot as plt

#Asteroid data
ast_x = 1879.6757
ast_y = 1526.6528

#Reference star data
#X, Y, RAh, RAm, RAs, DECdeg, DECmin, DECsec, mag
##ref_stars = [[58.594360, 14.510600, 20, 43, 18.647650, -5, 9, 52.895000, 9.12],
##[62.936470, 17.337510, 20, 42, 29.434250, -5, 18, 4.429000, 8.51],
##[65.792170, 29.402400, 20, 41, 58.070800, -5, 52, 27.624000, 9.96],
##[55.683940, 31.203150, 20, 43, 53.813550, -5, 57, 14.799000, 8.47],
##[53.172820, 28.121940, 20, 44, 22.179900, -5, 48, 24.180400, 9.29],
##[41.778100, 17.717280, 20, 46, 36.443850, -5, 18, 46.367000, 10.77]]

ref_stars = [[702.28, 488.0031, 264.611303, -25.0046098, 12.93],
             [2751.9015, 852.2215, 264.409387, -24.9691956, 12.423],
             [2052.0308, 895.0708, 264.478409, -24.9658256, 8.718],
             [1759.227, 1095.8338, 264.507497, -24.9494675, 9.977],
             [837.9692, 1452.1108, 264.599696, -24.918367, 11.844],
             [380.9108, 1716.3477, 264.645243, -24.8950537, 12.281]]

#Calculate b1, b2, a11, a21, a12, and a22
N = len(ref_stars)
star_ras = []
star_decs = []
for star in ref_stars:
##    star_ras.append(HMStoDeg(star[2], star[3], star[4]))
##    star_decs.append(DMStoDeg(star[5], star[6], star[7]))
    star_ras.append(star[2])
    star_decs.append(star[3])

sum_ras = sum(star_ras)
sum_decs = sum(star_decs)
sum_ra_x = 0
sum_ra_y = 0
sum_dec_x = 0
sum_dec_y = 0
for n in range(N):
    sum_ra_x += star_ras[n]*ref_stars[n][0]
    sum_ra_y += star_ras[n]*ref_stars[n][1]
    sum_dec_x += star_decs[n]*ref_stars[n][0]
    sum_dec_y += star_decs[n]*ref_stars[n][1]
    
sum_x = sum(get_col(ref_stars, 0))
sum_xsquares = sum([x**2 for x in get_col(ref_stars, 0)])
sum_y = sum(get_col(ref_stars, 1))
sum_ysquares = sum([x**2 for x in get_col(ref_stars, 1)])
sum_x_y = 0
for star in ref_stars:
    sum_x_y += star[0]*star[1]

matrix = [[N, sum_x, sum_y],
          [sum_x, sum_xsquares, sum_x_y],
          [sum_y, sum_x_y, sum_ysquares]]
ans_col_ra = [sum_ras, sum_ra_x, sum_ra_y]
ans_col_dec = [sum_decs, sum_dec_x, sum_dec_y]

##print(matrix)
##print(ans_col_ra)
##print(ans_col_dec)

ra_results = cramer(matrix, ans_col_ra)
dec_results = cramer(matrix, ans_col_dec)

b1 = ra_results[0]
a11 = ra_results[1]
a12 = ra_results[2]
b2 = dec_results[0]
a21 = dec_results[1]
a22 = dec_results[2]

# Define RA and DEC equations
def ra_eq(x, y):
    return b1 + a11*x + a12*y
def dec_eq(x, y):
    return b2 + a21*x + a22*y

# Calculate star coordinates
starRA = []
starDEC = []
for n in range(N):
    starRA.append(ra_eq(ref_stars[n][0], ref_stars[n][1]))
    starDEC.append(dec_eq(ref_stars[n][0], ref_stars[n][1]))

#Calculate residuals
starRA_residuals = []
starDEC_residuals = []
for n in range(N):
    starRA_residuals.append(star_ras[n] - ra_eq(ref_stars[n][0], ref_stars[n][1]))
    starDEC_residuals.append(star_decs[n] - dec_eq(ref_stars[n][0], ref_stars[n][1]))
##print("Star RA residuals:", starRA_residuals)
##print("Star DEC residuals:", starDEC_residuals)
##print()
#Calculate standard deviations
stdRA = sqrt(sum([x**2 for x in starRA_residuals]) / (N-3))
stdDEC = sqrt(sum([x**2 for x in starDEC_residuals]) / (N-3))

abs_std = sqrt(stdRA**2 + stdDEC**2)

#Calculate asteroid coordinates
astRA = RAdecimalToHMS(ra_eq(ast_x, ast_y))
print(stdRA)
stdRA = RAdecimalToHMS(stdRA)
astDEC = DECdecimalToDMS(dec_eq(ast_x, ast_y))
stdDEC = DECdecimalToDMS(stdDEC)

#Function for formatting hh:mm:ss
def format(matrix):
    return str(matrix[0]).zfill(2) + ":" + str(matrix[1]).zfill(2) + ":" + \
        str(matrix[2]).zfill(4)

print("RA: " + format(astRA) + " +/- " + format(stdRA))
print("DEC: " + format(astDEC) + " +/- " + format(stdDEC))
print()
##Visualization
##ast_size = abs_std*500
##size = [x[8]**3 for x in ref_stars]
##plt.scatter(starRA, starDEC, s=size)
##plt.plot(ra_eq(ast_x, ast_y), dec_eq(ast_x, ast_y), 'ro', ms = ast_size*2)
##plt.show()


