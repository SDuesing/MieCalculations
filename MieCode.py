# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 00:32:08 2019

@author: Sebastian Duesing
"""

import numpy as np
import PyMieScatt as ps

table = np.loadtxt(
    'C:/Users/duesing/Documents/PhD/MelCol/Optik/Calculations/Mie_R/20150626a_coarseMode_5micron_1particle_sd1.3.nsd',
    sep="\t", )
print(table)
diameter = table[0, :][4:]
concentrations = table[1::2, 4:]
numberOfScans = np.ma.shape(concentrations)[1]

wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
for i in range(0, numberOfScans - 1):
    concOfScan = concentrations[i, :][:]
    for w in wavelengths:
        bin = 0
        for d in diameter:
            result = ps.MieQ(1.53 + 0.015j, d, w, asDict=True)*np.pi*1.0/4.0*(d*1.0e-9.0)**2*concOfScan[bin]
            print(result)
            # Qback = result['Qback'] * konz / 1e-6 * 1e6 * pi * 1 / 4 * (diameter * 1e-9) ^ 2)
            # Qabs = result['Qabs']
            # Qsca = result['Qsca']
            # Qext = result['Qext']

            #np.trapz(y, x=diameter, axis=-1)

