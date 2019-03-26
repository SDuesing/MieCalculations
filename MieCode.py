# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 00:32:08 2019

@author: Sebastian Duesing
"""

import numpy as np
import PyMieScatt as ps
from integrate import integrateTrap
from hygroGrowth import wetdiameter


table = np.loadtxt(
    'test.txt',
    delimiter="\t")

diameter = table[0, :][4:]
concentrations: float = table[1::2, 4:]
numberOfScans = np.ma.shape(concentrations)[0]

wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
finalCoefficients = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4))
for i in range(0, numberOfScans, 1):
    concOfScan = concentrations[i, :][:]
    countW = 0
    for w in wavelengths:
        diameterBin = 0
        QbackVector = np.zeros((np.ma.shape(concOfScan)[0], 1))
        QextVector = np.zeros((np.ma.shape(concOfScan)[0], 1))

        QscaVector = np.zeros((np.ma.shape(concOfScan)[0], 1))
        QabsVector = np.zeros((np.ma.shape(concOfScan)[0], 1))
        for d in diameter:
            # print(d, w)
            numberConcentration = concOfScan[diameterBin]
            result = ps.MieQ(1.53 + 0.015j, d, w, asDict=True)
            QbackVector[diameterBin, 0] = result['Qback'] * np.pi * 1.0 / 4.0 * (
                    d * 1e-9) ** 2.0 * numberConcentration / (4.0 * np.pi)
            QabsVector[diameterBin, 0] = result['Qabs'] * np.pi * 1.0 / 4.0 * (d * 1e-9) ** 2.0 * numberConcentration
            QscaVector[diameterBin, 0] = result['Qsca'] * np.pi * 1.0 / 4.0 * (d * 1e-9) ** 2.0 * numberConcentration
            QextVector[diameterBin, 0] = result['Qext'] * np.pi * 1.0 / 4.0 * (d * 1e-9) ** 2.0 * numberConcentration

            diameterBin += 1
        SigmaAbs = integrateTrap(y=QabsVector[:, 0], x=np.ma.transpose(diameter))
        SigmaExt = integrateTrap(y=QextVector[:, 0], x=np.ma.transpose(diameter))
        SigmaSca = integrateTrap(y=QscaVector[:, 0], x=np.ma.transpose(diameter))
        SigmaBack = integrateTrap(y=QbackVector[:, 0], x=np.ma.transpose(diameter))

        columnStart = 0 + countW * 4
        columnEnd = 4 + countW * 4
        coefficients = np.array([[SigmaExt, SigmaSca, SigmaAbs, SigmaBack]])
        coefficients = coefficients/1e-6*1e6
        print(coefficients)
        finalCoefficients[i, columnStart:columnEnd] = coefficients
        countW += 1

np.savetxt(fname="results.txt", X=finalCoefficients, delimiter="\t", newline="\n")


kappa = np.linspace(0.1, 0.4, 40)
diameter = np.linspace(0.1, 2.0, 100)
result = np.zeros((diameter.shape[0]*kappa.shape[0], 3))

count = 0
for k in kappa:
    for d in diameter:
        d_wet = wetdiameter(d, k, 20, 0.6)
        row = np.array([k, d, d_wet])
        result[count, 0:3] = np.transpose(row)
        count += 1
print(result)
np.savetxt(fname='d_wet_table.txt', X=result, delimiter="\t", newline="\n")