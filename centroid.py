from math import *
import numpy as np
import matplotlib.pyplot as plt

def calc_centroid(matrix, x0, y0):
    x0 = 0
    y0 = 0
    count = 0
    weight_x = 0
    weight_y = 0
    for y in range(0, len(matrix)):
        row = matrix[y]
        for x in range(0, len(row)):
            count += row[x]
            weight_x += row[x]*(x - x0)
            weight_y += row[x]*(y - y0)
    xcm = weight_x / count
    ycm = weight_y / count
    print(xcm, ycm)
    return [[xcm, calc_uncertainty(matrix, xcm, ycm, x0, y0, count)[0]], [ycm, calc_uncertainty(matrix, xcm, ycm, x0, y0, count)[1]]]

def calc_uncertainty(matrix, xcm, ycm, x0, y0, count):
    uxsum = 0
    uysum = 0
    for y in range(0, len(matrix)):
        row = matrix[y]
        for x in range(0, len(row)):
            uxsum += ((((x - x0) - xcm) / count)**2 * row[x])
            uysum += ((((y - y0) - ycm) / count)**2 * row[x])
    delx = sqrt(uxsum)
    dely = sqrt(uysum)
    return [delx, dely]

# PART A
##pix_array = [[0, 33, 21, 33, 8], [0, 56, 51, 53, 26], [23, 120, 149, 73, 18], \
##             [55, 101, 116, 50, 16], [11, 78, 26, 2, 10]]
##calc_centroid(pix_array, 0, 0)


# PART B
##data2 = [[0,0,0,0,0,0,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0], [1,1,1,1,1,1,1], \
##         [0,0,0,0,0,0,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0]]
##np.savetxt('centroidMatrix.csv',data2, delimiter = ',')
##data = np.loadtxt('centroidMatrix.csv', dtype = float, delimiter=',')
##pix_array = data.tolist()
##calc_centroid(pix_array, 0, 0)


# PART C
pixels = [[10000,10000,10000], [20000,65000,20000], [6100,21000,13000]]
size = []
xcoords = []
ycoords = []
for y in range(0, len(pixels)):
    row = pixels[y]
    for x in range(0, len(row)):
        xcoords.append(x)
        ycoords.append(y)
        size.append(row[x]**0.9)

plt.scatter(xcoords, ycoords, s=size)
plt.plot(calc_centroid(pixels, 0, 0)[0][0], calc_centroid(pixels, 0, 0)[1][0], 'rX')
plt.show()
