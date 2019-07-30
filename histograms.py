import numpy as np
import matplotlib.pyplot as plt
from odlib import *

data = np.loadtxt('elements.csv', dtype=float, delimiter=',')
elements = data.tolist()

el = elements[0]
avg = mean(el)
stde = stdev(el)
plt.hist(el, bins = np.arange(avg-stde*3, avg+stde*3, stde*0.5), color='white', edgecolor='black')
plt.title("Semimajor axis")
plt.ylabel("count")
plt.show()
