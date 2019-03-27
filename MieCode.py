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
import random

# setup
mixture = "CS"  # CS for core-shell
wet = True
if wet:
    print("wet")
else:
    print("dry")
filePathSizeDist = 'X:/p3_working/SDuesing-home/Data/MelCol/20150628/20150628b_lastcorr_smoothed_g0_diffkorr_ab8nm_totKorr_cpceff_inletcorr_final.nsd'
filePathActos = 'X:\\p3_working\\SDuesing-home\\Data\\MelCol\\20150628\\20150628b_ver2.txt'

fVolBC = 0.05  # 5% black carbon content
numberConcentrationError: float = 0.1  # 10 percent error size distribution

monteCarloIterations = 50

# refractive indices
realPartWater = 1.33
realPartBC = 1.75
realPartSol = 1.53

imPartWater = 1e-6
imPartBc = 0.55
imPartSol = 1e-6

RealPartSolError = 0.5  # percent
RealPartBCError = 4.  # percent
RealPartWaterError = 0.5  # percent

imPartSolError = 0.  # percent
imPartBCError = 6.6  # percent
imPartWaterError = 0.  # percent

# setup end


tableSizeDist = np.genfromtxt(
    filePathSizeDist,
    delimiter="\t")

file = os.path.basename(filePathSizeDist)
prefix = file.split(sep=".")[0]

diameter = tableSizeDist[0, :][4:]
concentrations: float = tableSizeDist[1::2, 4:]
ScanTimes = tableSizeDist[1::2, 0]
numberOfScans = np.ma.shape(concentrations)[0]
wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
finalCoefficients = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4), dtype=float)
finalCoefficientsSd = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4), dtype=float)
my_list = list(range(100))

numberOfCalculations = np.ma.shape(concentrations)[0] * diameter.size * wavelengths.size * monteCarloIterations
print(numberOfCalculations)
n = 0.

for i in range(0, numberOfScans, 1):
    concOfScan = concentrations[i, :][:]
    SigmaBackVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
    SigmaExtVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
    SigmaScaVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
    SigmaAbsVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
    for monte in range(monteCarloIterations):
        conc = concOfScan * (1. + random.uniform(-1., 1.) * numberConcentrationError)
        # generate normal distributed RH from measured values
        refIndBC = np.random.normal(realPartBC, realPartBC * RealPartBCError / 100., 1) + np.random.normal(
            imPartBc, imPartBc * imPartBCError / 100., 1) * 1j
        refIndSol = np.random.normal(realPartSol, realPartSol * RealPartSolError / 100., 1) + np.random.normal(
            imPartSol, imPartSol * imPartSolError / 100., 1) * 1j
        if wet:
            refIndWater = np.random.normal(realPartWater, realPartWater * RealPartWaterError / 100.,
                                           1) + np.random.normal(
                imPartWater, imPartWater * imPartWaterError / 100., 1) * 1j

            # generate normal distributed RH from measured values
            meanRH = 0.6
            sdRH = 0.05
            relHum = np.random.normal(meanRH, sdRH, 1)

            # generate normal distributed T from measured values
            meanT = 20
            sdT = 5
            T = np.random.normal(meanT, sdT, 1)


        #print(conc[0])



        QBackVector = np.zeros((np.ma.shape(conc)[0], wavelengths.size), dtype=float)
        QExtVector = np.zeros((np.ma.shape(conc)[0], wavelengths.size), dtype=float)
        QScaVector = np.zeros((np.ma.shape(conc)[0], wavelengths.size), dtype=float)
        QAbsVector = np.zeros((np.ma.shape(conc)[0], wavelengths.size), dtype=float)

        for countW in range(wavelengths.size):
            diameterBin = 0
            for d in diameter:
                # print(d, w)
                numberConcentration: float = conc[diameterBin]
                diameterBC = d * fVolBC ** (1. / 3.)
                if wet:
                    d_wet = wetdiameter(kappa=0.3, rh=relHum, ddry=d, temperature=T)
                    dParticle = d_wet
                    volumeWater = 1. / 6. * np.pi * (math.pow(d_wet, 3.) - math.pow(d, 3.))
                    volumeSol = np.pi * 1. / 6. * (math.pow(d, 3.) - math.pow(diameterBC, 3.))
                    volumeBC = fVolBC * 1. / 6. * np.pi * math.pow(d, 3.)
                    totalVolume = np.pi * 1. / 6. * math.pow(d_wet, 3.)
                    # print("Total Volume - Sum of volumes: ", totalVolume - volumeWater - volumeSol - volumeBC)
                    if mixture == "CS":
                        fout = prefix + "_" + str(monteCarloIterations) + "_wet_CS"
                        mShell = volumeSol / (volumeWater + volumeSol) * refIndSol + (
                                1 - volumeSol / (volumeWater + volumeSol)) * refIndWater
                        mCore = refIndBC
                        # print(mShell)
                    elif mixture == "HOM":
                        fout = prefix + "_" + str(monteCarloIterations) + "_wet_HOM"
                        fVoldry = wetVolumefraction(kappa=0.3, rh=0.7, ddry=d, temperature=15)
                        refIndAerosol = fVolBC * refIndBC + (1. - fVolBC) * refIndSol
                        refIndWet: Union[complex, Any] = refIndWater * (1. - fVoldry) + refIndAerosol * fVoldry
                        refInd = refIndWet

                else:
                    dParticle = d
                    if mixture == "CS":
                        fout = prefix + "_" + str(monteCarloIterations) + "_dry_CS"
                        mShell = refIndSol
                        mCore = refIndBC

                    elif mixture == "HOM":
                        refInd = fVolBC * refIndBC + (1. - fVolBC) * refIndSol
                        fout = prefix + "_" + str(monteCarloIterations) + "_dry_HOM"

                if mixture == "CS":
                    result = ps.MieQCoreShell(dCore=diameterBC, mShell=mShell, mCore=mCore, dShell=d,
                                              wavelength=wavelengths[countW], asDict=True)
                elif mixture == "HOM":
                    result = ps.MieQ(refInd, diameter=dParticle, wavelength=wavelengths[countW], asDict=True)

                QBackVector[diameterBin, countW] = result['Qback'] * np.pi * 1.0 / 4.0 * (
                        dParticle * 1e-9) ** 2.0 * numberConcentration / (4.0 * np.pi)
                QAbsVector[diameterBin, countW] = result['Qabs'] * np.pi * 1.0 / 4.0 * (
                        dParticle * 1e-9) ** 2.0 * numberConcentration
                QScaVector[diameterBin, countW] = result['Qsca'] * np.pi * 1.0 / 4.0 * (
                        dParticle * 1e-9) ** 2.0 * numberConcentration
                QExtVector[diameterBin, countW] = result['Qext'] * np.pi * 1.0 / 4.0 * (
                        dParticle * 1e-9) ** 2.0 * numberConcentration
                #print(QExtVector)
                drawProgressBar(n / numberOfCalculations, barLen=20)
                n += 1
                diameterBin += 1
                #countW += 1

            for nWavelengths in range(wavelengths.size):
                SigmaAbs = integrateTrap(y=QAbsVector[:, nWavelengths], x=np.ma.transpose(diameter))
                SigmaExt = integrateTrap(y=QExtVector[:, nWavelengths], x=np.ma.transpose(diameter))
                SigmaSca = integrateTrap(y=QScaVector[:, nWavelengths], x=np.ma.transpose(diameter))
                SigmaBack = integrateTrap(y=QBackVector[:, nWavelengths], x=np.ma.transpose(diameter))
                SigmaBackVectorMonte[monte, nWavelengths] = SigmaBack
                SigmaExtVectorMonte[monte, nWavelengths] = SigmaExt
                #print("SigmaExtMonteMatrix Monte durchlauf" + str(monte))
                #print(SigmaExtVectorMonte)
                SigmaScaVectorMonte[monte, nWavelengths] = SigmaSca
                SigmaAbsVectorMonte[monte, nWavelengths] = SigmaAbs

    coefficients = []
    coefficientsSd = []
    for nWavelengths in range(wavelengths.size):
        SigmaExtMean = np.mean(SigmaExtVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficients.append(SigmaExtMean)
        SigmaScaMean = np.mean(SigmaScaVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficients.append(SigmaScaMean)
        SigmaAbsMean = np.mean(SigmaAbsVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficients.append(SigmaAbsMean)
        SigmaBackMean = np.mean(SigmaBackVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficients.append(SigmaBackMean)

        SigmaExtMean = np.std(SigmaExtVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficientsSd.append(SigmaExtMean)
        SigmaScaMean = np.std(SigmaScaVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficientsSd.append(SigmaScaMean)
        SigmaAbsMean = np.std(SigmaAbsVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficientsSd.append(SigmaAbsMean)
        SigmaBackMean = np.std(SigmaBackVectorMonte[:, nWavelengths]) / 1.e-6 * 1.e6
        coefficientsSd.append(SigmaBackMean)
    #print(coefficients)
    finalCoefficients[i, :] = coefficients
    finalCoefficientsSd[i, :] = coefficients


np.savetxt(fname=fout+"_SD.txt", X=finalCoefficientsSd, delimiter="\t")
np.savetxt(fname=fout+".txt", X=finalCoefficients, delimiter="\t")

plt.plot(ScanTimes, finalCoefficients[:, 0], "bo")
plt.show()
