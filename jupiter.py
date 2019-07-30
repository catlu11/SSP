import numpy as np
from math import *
# Jupiter data rotation

j_x = 1784.468272
j_y = 1228.092276
theta = -33

x = np.array([1269.488953, 2023.832215, 2634.376353])
y = np.array([1611.401073, 1033.836061, 603.464293])

# Translate coordinates
x -= j_x
y -= j_y

# Determine x' 
xprime = [cos(-theta*pi/180)*x[i] + sin(-theta*pi/180)*y[i] for i in range(len(x))]

for i in xprime:
    print(i)
