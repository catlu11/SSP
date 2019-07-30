from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from od6 import *

plt.figure()
ax = plt.axes(projection='3d')

# Plot sun
ax.scatter3D(0,0,0, color="#F39C12")

# Generate times
earth_times = np.arange(2458683.774, 2459300.774, 1)
ast_times = np.arange(2458683.774, 2460400.774, 1)
mars_times = np.arange(2458683.774, 2460000.774, 1)

# Plot earth (ecliptic coordinates)
earth_elements = [1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 189.5127]
e = []
for t in earth_times:
    e.append(generate_coords(earth_elements, t))
e_x = [x[0] for x in e]
e_y = [x[1] for x in e]
e_z = [x[2] for x in e]
ax.scatter3D(e_x, e_y, e_z, s=1, c='g')

# Plot asteroid (ecliptic coordinates)
ast_elements = [2.64472, 0.40668, 9.53748, 280.97715, 356.76770, 1.74872]
a = []
for t in ast_times:
    a.append(generate_coords(ast_elements, t))
a_x = [x[0] for x in a]
a_y = [x[1] for x in a]
a_z = [x[2] for x in a]
ax.scatter3D(a_x, a_y, a_z, s=15, c='b')

# Plot Mars (ecliptic coordinates)
mars_elements = [1.524, 0.0934, 1.85061, 49.57854, 336.04084, 158.100]
m = []
for t in mars_times:
    m.append(generate_coords(mars_elements, t))
m_x = [x[0] for x in m]
m_y = [x[1] for x in m]
m_z = [x[2] for x in m]
ax.scatter3D(m_x, m_y, m_z, s=1, c='r')

plt.title("Orbit of 2002 KM6 around the Sun")
ax.set_xlabel("x (AU)")
ax.set_ylabel("y (AU)")
ax.set_zlabel("z (AU)")
plt.show()
