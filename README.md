# SSP
A collection of Python programs that I wrote for my orbit determination research at the Summer Science Program.
* odlib.py - library containing functions for vector operations and celestial coordinate conversions
* od5.py - calculates the orbital elements of any object rotating around the Sun, given its right ascension and declination on three separate days (time in Julian days) and the equatorial Earth-Sun vectors on those days (in AU)
* od6.py - generates the Sun-object vector at any time given the object’s orbital elements 
* lspr.py - calculates the coordinates of an object in an image given the coordinates of reference bodies visible around the object
* centroid.py -  contains functions for calculating the centroid of a matrix
* histograms.py - generates a histogram for a list of values
* monte_carlo.py - runs a Monte Carlo simulation to sample possible orbital element values, given their error ranges
* residuals.py - calculates the residuals of observed and calculated coordinates
* radec_to_altazi.py - converts right ascension and declination coordinates to altitude and azimuth values, given the time (LST) when the coordinates were taken
* 3dplot.py - generates a 3D plot of an object’s orbit around the Sun given the object's orbital elements 
