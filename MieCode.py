# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 00:32:08 2019

@author: Sebastian Duesing
"""
from typing import Any, Union
import PyMieScatt as ps
from integrate import integrateTrap
from hygroGrowth import *
import matplotlib.pyplot as plt

# setup
mixture = "HOM" # CS for core-shell
refIndAerosol = 1.53 + 0.015j
refIndWater = 1.33 + 1e-6j
wet = True
if wet:
    print("wet")
else:
    print("dry")
# setup end
table = np.genfromtxt(
    'test.txt',
    delimiter="\t")
diameter = table[0, :][4:]
concentrations: float = table[1::2, 4:]
numberOfScans = np.ma.shape(concentrations)[0]
wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
finalCoefficients = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4), dtype=float)

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
            numberConcentration: float = concOfScan[diameterBin]
            if wet:
                # print("wet calculation")
                fout = "results_wet.txt"
                fVolWater = wetVolumefraction(kappa=0.3, rh=0.7, ddry=d, temperature=15)
                refIndWet: Union[complex, Any] = refIndWater * (1. - fVolWater) + refIndAerosol * fVolWater
                d_wet = wetdiameter(kappa=0.3, rh=0.7, ddry=d, temperature=15)
                d = d_wet
                refInd = refIndWet
            else:
                refInd = refIndAerosol
                fout = "results_dry.txt"
            if mixture=="CS":
                # result = ps.MieQCoreShell(dCore=d*fVolBC**(1./3.)refInd, d, w, asDict=True)
            elif mixture == "HOM":
                result = ps.MieQ(refInd, d, w, asDict=True)
            QbackVector[diameterBin, 0] = result['Qback'] * np.pi * 1.0 / 4.0 * (
                    d * 1e-9) ** 2.0 * numberConcentration / (4.0 * np.pi)
            QabsVector[diameterBin, 0] = result['Qabs'] * np.pi * 1.0 / 4.0 * (
                    d * 1e-9) ** 2.0 * numberConcentration
            QscaVector[diameterBin, 0] = result['Qsca'] * np.pi * 1.0 / 4.0 * (
                    d * 1e-9) ** 2.0 * numberConcentration
            QextVector[diameterBin, 0] = result['Qext'] * np.pi * 1.0 / 4.0 * (
                    d * 1e-9) ** 2.0 * numberConcentration

            diameterBin += 1
        SigmaAbs = integrateTrap(y=QabsVector[:, 0], x=np.ma.transpose(diameter))
        SigmaExt = integrateTrap(y=QextVector[:, 0], x=np.ma.transpose(diameter))
        SigmaSca = integrateTrap(y=QscaVector[:, 0], x=np.ma.transpose(diameter))
        SigmaBack = integrateTrap(y=QbackVector[:, 0], x=np.ma.transpose(diameter))

        columnStart = 0 + countW * 4
        columnEnd = 4 + countW * 4
        coefficients = np.array([[SigmaExt, SigmaSca, SigmaAbs, SigmaBack]])
        coefficients = coefficients / 1e-6 * 1e6

        finalCoefficients[i, columnStart:columnEnd] = coefficients
        countW += 1

    np.savetxt(fname=fout, X=finalCoefficients, delimiter="\t")

kappa = np.linspace(0.1, 0.4, 40)
diameter = np.linspace(0.1, 5.0, 500)
result = np.zeros((diameter.shape[0] * kappa.shape[0], 3))

count = 0
for k in kappa:
    for d in diameter:
        d_wet = wetdiameter(d, k, 20, 0.9) / d
        row = np.array([k, d, d_wet])
        result[count, 0:3] = np.transpose(row)
        count += 1
np.savetxt(fname='d_wet_table.txt', X=result, delimiter="\t")

x = result[:, 0]
y = result[:, 1]
z = result[:, 2]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cmhot = plt.get_cmap("hot")
c = z
ax.scatter(x, y, z, zdir='z', c=c, cmap=cmhot)
plt.show()
