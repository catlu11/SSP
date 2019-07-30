from odlib import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

RA9 = [17.66471, 17.6647348, 17.66517685, 17.66633422, 17.66957495, 17.66989513, \
       17.67145873, 17.67266549, 17.67395345, 17.67591804, 17.67681085, 17.66377112, \
       17.66501181, 17.66525481, 17.6657882, 17.66616314, 17.66879573, 17.66907753, \
       17.66965404, 17.67243804]
RAJ = [17.664722, 17.664743, 17.665176, 17.666338, 17.669568, 17.669895, 17.671452, \
       17.67267, 17.673959, 17.675919, 17.676814, 17.663778, 17.665017, 17.665261, \
       17.665794, 17.666175, 17.66881, 17.669083, 17.669663, 17.672449]

DEC9 = [-27.8328625, -27.8441506, -27.8750767, -27.8079114, -27.9252837, -27.8687625,\
        -27.9063159, -27.8057003, -27.8071242, -27.8636878, -27.8393767, -27.7806937,\
        -27.6811923, -27.7716092, -27.7873273, -27.6734623, -27.6803645, -27.7971312,\
        -27.7436503, -27.7564373]
DECJ = [-27.832693, -27.843905, -27.874865, -27.807628, -27.925135, -27.868543, -27.906242,\
        -27.805697, -27.807262, -27.863801, -27.839808, -27.780385, -27.681199, -27.771325,\
        -27.787152, -27.673465, -27.680434, -27.796966, -27.743515,-27.756456]
starRA_residuals = []
starDEC_residuals = []
for n in range(len(RA9)):
    starRA_residuals.append(RAJ[n] - RA9[n])
    starDEC_residuals.append(DECJ[n] - DEC9[n])
starRA_residuals = np.array(starRA_residuals)
starDEC_residuals = np.array(starDEC_residuals)

avgRA = mean(starRA_residuals)*3600*15
stdeRA = stdev(starRA_residuals)*3600*15
stdmRA = stdeRA/sqrt(len(RA9))
avgDEC = mean(starDEC_residuals)*3600
stdeDEC = stdev(starDEC_residuals)*3600
stdmDEC = stdeDEC/sqrt(len(RA9))

##plt.hist(starRA_residuals*3600*15, bins = np.arange(avgRA-stdeRA*3, avgRA+stdeRA*3, stdeRA*0.5), color='white', edgecolor='black')
##plt.title("RA residuals")
##plt.xlabel("deviation (arcseconds)")
plt.hist(starDEC_residuals*3600, bins = np.arange(avgDEC-stdeDEC*3, avgDEC+stdeDEC*3, stdeDEC*0.5), color='white', edgecolor='black')
plt.title("DEC residuals")
plt.xlabel("deviation (arcseconds)")
plt.show()
          
print("RA Mean:", avgRA)
print("RA Standard Deviation:", stdeRA)
print("RA Standard Deviation of Mean:", stdmRA)
print()
print("DEC Mean:", avgDEC)
print("DEC Standard Deviation:", stdeDEC)
print("DEC Standard Deviation of Mean:", stdmDEC)
