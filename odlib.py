from math import *
from copy import deepcopy

# Calculates the dot product of two vectors: vector1 and vector2.
def dot(vector1, vector2):
    min_len = min(len(vector1), len(vector2))
    dot_p = 0
    for i in range(0, min_len):
        dot_p += (vector1[i] * vector2[i])
    return dot_p


# Calculates the cross product of two vectors: vector 1 and vector2.
def cross(vector1, vector2):
    x_coef = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1])
    y_coef = (vector1[0] * vector2[2]) - (vector1[2] * vector2[0])
    z_coef = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])
    return [x_coef, -y_coef, z_coef]


# Calculates the triple product of three vectors: vector1, vector2, and vector3.
def tri_product(vector1, vector2, vector3):
    cross23 = cross(vector2, vector3)
    return dot(vector1, cross23)


# Calculates mean of a list of numbers.
def mean(num_list):
    return sum(num_list) / len(num_list)


# Calculates standard deviation of a list of numbers.
def stdev(num_list):
    avg = mean(num_list)
    terms = [(x-avg)**2 for x in num_list]
    temp_sum = sum(terms)
    return sqrt(temp_sum / (len(num_list)-1))


# Converts RA in HMS to an angle in decimalized degrees.
def HMStoDeg(hours, minutes, seconds):
    return hours*15 + minutes*(0.25) + seconds*(0.25/60)


# Converts declination in DMS to an angle in degrees.
def DMStoDeg(degrees, minutes, seconds):
    ans = abs(degrees) + minutes / 60 + seconds / 3600
    if(degrees > 0):
        return ans
    return -ans


# Displays the RA in HMS, given the decimalized value in degrees.
def RAdecimalToHMS(angle):
    hours = int(angle / 15)
    angle %= 15
    minutes = int(angle / 0.25)
    angle %= 0.25
    seconds = angle / (0.25/60)
    return [hours, minutes, seconds]


# Displays the declination in DMS given the decimalized value in degrees.
def DECdecimalToDMS(angle):
    degrees = int(angle)
    angle -= degrees
    if(angle < 0):
        angle *= -1
    minutes = int(angle / (1/60))
    angle %= (1/60)
    seconds = angle / (1/3600)
    return [degrees, minutes, seconds]
    

# Calculates the magnitude of a vector.
def mag(vector):
    return sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)


# Given the sine and cosine of an angle, returns the value of the angle in radians.
def trig_to_rad(sine, cosine):
    if(sine < 0):
        return 2*pi - acos(cosine)
    elif(cosine < 0):
        return pi - asin(sine)
    else:
        return asin(sine)


# Rotates a matrix theta degrees CCW about the y axis, phi degrees CW about x, and epsi degrees CW about z.
def rotation_sequence(vector, theta, phi, epsi):
    vprime = vector
    vprime = rotate_z(rotate_x(rotate_y(vprime, theta), -phi), -epsi)
    return vprime


# Rotates a matrix theta degrees CCW about the y-axis.
def rotate_y(vector, theta):
    theta *= pi/180
    x = vector[0]*cos(theta) - vector[2]*sin(theta)
    z = vector[2]*cos(theta) + vector[0]*sin(theta)
    return [x, vector[1], z]


# Rotates a matrix phi degrees CCW about the x-axis.
def rotate_x(vector, phi):
    phi *= pi/180
    y = vector[1]*cos(phi) + vector[2]*sin(phi)
    z = vector[2]*cos(phi) - vector[1]*sin(phi)
    return [vector[0], y, z]


# Rotates a matrix epsi degrees CCW about the z-axis.
def rotate_z(vector, epsi):
    epsi *= pi/180
    x = vector[0]*cos(epsi) - vector[1]*sin(epsi)
    y = vector[1]*cos(epsi) + vector[0]*sin(epsi)
    return [x, y, vector[2]]


# Returns the determinant of a 2x2 matrix
def det2x2(matrix):
    return (matrix[0][0]*matrix[1][1]) - (matrix[0][1]*matrix[1][0])


# Returns the determinant of a 3x3 matrix
def det3x3(matrix):
    dets = []
    for i in range(len(matrix)):
        coef = matrix[0][i]
        row1 = []
        row2 = []
        for j in [x for x in range(len(matrix)) if x != i]:
            row1.append(matrix[1][j])
            row2.append(matrix[2][j])
        two_by_two = [row1, row2]
        dets.append(coef*det2x2(two_by_two))
    return dets[0] - dets[1] + dets[2]


# Replaces the column of a matrix with a list of coefficients
def spliceColumn(matrix, coeff, col):
    n_matrix = deepcopy(matrix)
    for r in range(len(matrix)):
        n_matrix[r][col] = coeff[r]
    return n_matrix


# Calculates the solution of three equations using Cramer's rule
def cramer(matrix, constants):
    d = det3x3(matrix)
    dx = det3x3(spliceColumn(matrix, constants, 0))
    dy = det3x3(spliceColumn(matrix, constants, 1))
    dz = det3x3(spliceColumn(matrix, constants, 2))
    return [dx/d, dy/d, dz/d]


# Returns the column of an array in the form of a list
def get_col(array, x):
    col = []
    for row in array:
        col.append(row[x])
    return col


# Returns a list of percent errors for corresponding estimate and prediction lists
def calc_percent_error(expected, actual):
    errors = []
    for i in range(len(expected)):
        errors.append(abs((actual[i] - expected[i])/expected[i]))
    return errors


















