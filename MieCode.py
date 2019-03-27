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
import math
import os
from progressbar import drawProgressBar

# setup
mixture = "CS"  # CS for core-shell
wet = True

refIndBC = 1.75 + 0.55j
refIndSol = 1.53 + 1e-6j
refIndWater = 1.33 + 1e-6j

fVolBC = 0.05  # 5% black carbon content

if wet:
    print("wet")
else:
    print("dry")
# setup end
filePathSizeDist = 'X:/p3_working/SDuesing-home/Data/MelCol/20150628/20150628b_lastcorr_smoothed_g0_diffkorr_ab8nm_totKorr_cpceff_inletcorr_final.nsd'
filePathActos = 'X:\\p3_working\\SDuesing-home\\Data\\MelCol\\20150628\\20150628b_ver2.txt'

tableSizeDist = np.genfromtxt(
    filePathSizeDist,
    delimiter="\t")

file = os.path.basename(filePathSizeDist)
prefix = file.split(sep=".")[0]

diameter = tableSizeDist[0, :][4:]
concentrations: float = tableSizeDist[1::2, 4:]
numberOfScans = np.ma.shape(concentrations)[0]
wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
finalCoefficients = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4), dtype=float)
my_list = list(range(100))

numberOfCalculations = np.ma.shape(concentrations)[0] * diameter.size * wavelengths.size
print(numberOfCalculations)
n = 0.
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
            diameterBC = d * fVolBC ** (1. / 3.)
            if wet:
                d_wet = wetdiameter(kappa=0.3, rh=0.7, ddry=d, temperature=15)
                dParticle = d_wet
                volumeWater = 1. / 6. * np.pi * (math.pow(d_wet, 3.) - math.pow(d, 3.))
                volumeSol = np.pi * 1. / 6. * (math.pow(d, 3.) - math.pow(diameterBC, 3.))
                volumeBC = fVolBC * 1. / 6. * np.pi * math.pow(d, 3.)
                totalVolume = np.pi * 1. / 6. * math.pow(d_wet, 3.)
                # print("Total Volume - Sum of volumes: ", totalVolume - volumeWater - volumeSol - volumeBC)
                if mixture == "CS":
                    fout = prefix + "_wet_CS.txt"
                    mShell = volumeSol / (volumeWater + volumeSol) * refIndSol + (
                            1 - volumeSol / (volumeWater + volumeSol)) * refIndWater
                    mCore = refIndBC
                    # print(mShell)
                elif mixture == "HOM":
                    fout = prefix + "_wet_HOM.txt"
                    fVoldry = wetVolumefraction(kappa=0.3, rh=0.7, ddry=d, temperature=15)
                    refIndAerosol = fVolBC * refIndBC + (1. - fVolBC) * refIndSol
                    refIndWet: Union[complex, Any] = refIndWater * (1. - fVoldry) + refIndAerosol * fVoldry
                    refInd = refIndWet

            else:
                dParticle = d
                if mixture == "CS":
                    fout = prefix + "_dry_CS.txt"
                    mShell = refIndSol
                    mCore = refIndBC

                elif mixture == "HOM":
                    refInd = fVolBC * refIndBC + (1. - fVolBC) * refIndSol
                    fout = prefix + "_dry_HOM.txt"

            if mixture == "CS":
                result = ps.MieQCoreShell(dCore=diameterBC, mShell=mShell, mCore=mCore, dShell=d,
                                          wavelength=w, asDict=True)
            elif mixture == "HOM":
                result = ps.MieQ(refInd, diameter=dParticle, wavelength=w, asDict=True)

            QbackVector[diameterBin, 0] = result['Qback'] * np.pi * 1.0 / 4.0 * (
                    dParticle * 1e-9) ** 2.0 * numberConcentration / (4.0 * np.pi)
            QabsVector[diameterBin, 0] = result['Qabs'] * np.pi * 1.0 / 4.0 * (
                    dParticle * 1e-9) ** 2.0 * numberConcentration
            QscaVector[diameterBin, 0] = result['Qsca'] * np.pi * 1.0 / 4.0 * (
                    dParticle * 1e-9) ** 2.0 * numberConcentration
            QextVector[diameterBin, 0] = result['Qext'] * np.pi * 1.0 / 4.0 * (
                    dParticle * 1e-9) ** 2.0 * numberConcentration
            drawProgressBar(n/numberOfCalculations, barLen=20)
            n += 1
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

plt.plot(finalCoefficients[:, 0], "bo")
plt.show()
