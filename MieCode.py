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
from extractAmbientConditions import extractAmbient
import datetime as dt
import time
import pandas as pd
from averageKappaCyrielle import averageKappa, findKappa
import glob

# setup
mixtures = ["CS", "HOM"]
states = ["wet", "dry"]

fileChemieLaurent = "tabelle_masse_volfrac_refind_kappa.txt"
fCyrielle = "kappa_cyrielle.txt"

fList = glob.glob("C:\\Users\\duesing\\PycharmProjects\\MieCalculations\\ACTOSNSD\\*.nsd")
fListActos = glob.glob("C:\\Users\\duesing\\PycharmProjects\\MieCalculations\\ACTOSFiles\\*.txt")
# fList = glob.glob("*.nsd")
# fListActos = glob.glob("*.dat")

chemManual = False
kappaManual = False
fVolBC = 0.012  # 1.2% black carbon content
numberConcentrationError: float = 0.1  # 10 percent error size distribution

monteCarloIterations = 50

# refractive indices
realPartWater = 1.33
realPartBC = 1.75
realPartSol = 1.53

imPartWater = 1e-6
imPartBc = 0.55
imPartSol = 1e-6
# deviations
RealPartSolError = 0.5  # percent
RealPartBCError = 4.  # percent
RealPartWaterError = 0.5  # percent

imPartSolError = 0.  # percent
imPartBCError = 6.6  # percent
imPartWaterError = 0.  # percent

scanDuration = 122.  # time of one scan

colRH = 11
colTemp = 7
colTime = 1
colVolBC = 13

dHTDMA = [35, 50, 75, 110, 165, 265]

# setup end
# CS for core-shell
for m in range(np.size(mixtures)):
    mixture = mixtures[m]
    for s in range(np.size(states)):
        state = states[s]
        if state == "wet":
            print("wet")
        else:
            print("dry")

        for f in range(np.size(fList)):
            start = time.time()
            file = os.path.basename(fList[f])
            prefix = file.split(sep=".")[0]
            tableSizeDist = np.genfromtxt(
                fList[f],
                delimiter="\t")
            date = file[0:8]
            datum = dt.datetime.strptime(date, "%Y%m%d")
            actosData = np.genfromtxt(fListActos[f], delimiter='\t')
            while (actosData[0, colTime] > tableSizeDist[1, 0]) | (actosData[-1, colTime] < tableSizeDist[-1, 0]):
                if actosData[0, colTime] > tableSizeDist[1, 0]:
                    tableSizeDist = tableSizeDist[2:, :]
                if actosData[-1, colTime] < tableSizeDist[-1, 0]:
                    tableSizeDist = tableSizeDist[:-2, :]
            datesActos = [dt.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])) + dt.timedelta(seconds=x) for x in
                          actosData[:, colTime]]
            #print(dt.datetime.date(datesActos))
            diameter = tableSizeDist[0, :][4:]
            concentrations: float = tableSizeDist[1::2, 4:]
            ScanTimes = tableSizeDist[1::2, 0]
            datesScans = [dt.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8])) + dt.timedelta(seconds=x) for x in
                          ScanTimes]
            numberOfScans = np.ma.shape(concentrations)[0]
            print("Number of Scans: " + str(numberOfScans))
            wavelengths = np.array([355, 450, 525, 532, 624, 630, 880, 1064])
            finalCoefficients = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4 + 1), dtype=float)
            finalCoefficientsSd = np.zeros((numberOfScans, np.ma.shape(wavelengths)[0] * 4 + 1), dtype=float)
            my_list = list(range(100))

            numberOfCalculations = np.ma.shape(concentrations)[
                                       0] * diameter.size * wavelengths.size * monteCarloIterations
            print(numberOfCalculations)
            n = 0.

            if not chemManual:
                tableChemistryLaurent = pd.read_csv(
                    fileChemieLaurent,
                    sep="\t", header=None)
                datesChem = tableChemistryLaurent[0]
                dateChemistry = [dt.datetime.strptime(date, "%d.%m.%Y %H:%M") for date in datesChem]
                chemistryStart = int(np.argwhere(np.array(dateChemistry) > datesActos[0])[0])
                print(str(datesActos[-1]))
                chemistryEnd = int(np.argwhere(np.array(dateChemistry) <= datesActos[-1])[-1])
                volBCMean = np.mean(tableChemistryLaurent[colVolBC - 1][chemistryStart:(chemistryEnd + 1)])
                sdvolBC = np.std(tableChemistryLaurent[colVolBC - 1][chemistryStart:(chemistryEnd + 1)])
                print("volBC:" + str(volBCMean) + " SD: " + str(sdvolBC))


            if (not kappaManual) and state == "wet":
                fileCyrielle = fCyrielle
                dHtdma = np.array([35, 50, 75, 110, 165, 265])
                tableKappaCyr = pd.read_csv(
                    fileCyrielle,
                    sep="\t")
                kappaArray = averageKappa(tableKappaCyr, str(date))

            for i in range(0, numberOfScans, 1):
                timeOfScan = ScanTimes[i]
                concOfScan = concentrations[i, :][:]
                SigmaBackVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
                SigmaExtVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
                SigmaScaVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
                SigmaAbsVectorMonte = np.zeros((monteCarloIterations, wavelengths.size), dtype=float)
                for monte in range(monteCarloIterations):
                    volBC = np.random.uniform(volBCMean-sdvolBC, volBCMean+sdvolBC, 1)
                    if volBC < 0:
                        print(str(volBC))
                    conc = concOfScan * (1. + random.uniform(-1., 1.) * numberConcentrationError)
                    # generate normal distributed RH from measured values
                    refIndBC = np.random.normal(realPartBC, realPartBC * RealPartBCError / 100., 1) + np.random.normal(
                        imPartBc, imPartBc * imPartBCError / 100., 1) * 1j
                    refIndSol = np.random.normal(realPartSol, realPartSol * RealPartSolError / 100.,
                                                 1) + np.random.normal(
                        imPartSol, imPartSol * imPartSolError / 100., 1) * 1j
                    if state == "wet":
                        meanRH, sdRH, meanT, sdT = extractAmbient(actosData, rhCol=colRH, timeCol=colTime,
                                                                  tempCol=colTemp,
                                                                  startTime=ScanTimes[i],
                                                                  duration=scanDuration)
                        refIndWater = np.random.normal(realPartWater, realPartWater * RealPartWaterError / 100.,
                                                       1) + np.random.normal(
                            imPartWater, imPartWater * imPartWaterError / 100., 1) * 1j

                        # generate normal distributed RH from measured values
                        relHum = np.random.uniform(meanRH / 100. - sdRH / 100., meanRH / 100. + sdRH / 100., 1)

                        # generate normal distributed T from measured values
                        T = np.random.uniform(meanT - sdT, meanT + sdT, 1)

                        # generate normal distributed kappa from measured values
                        if kappaManual:
                            meanKappa = 0.366
                            sdKappa = 0.0121
                            Kappa = np.random.normal(meanKappa, sdKappa, 1)

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
                            if state == "wet":
                                if not kappaManual:
                                    meanKappa = findKappa(d, kappas=kappaArray[0:6], diameter=dHtdma)
                                    sdKappa = findKappa(d, kappas=kappaArray[6:], diameter=dHtdma)
                                    # Kappa = np.random.normal(meanKappa, sdKappa, 1)
                                    Kappa = np.random.uniform(meanKappa - sdKappa, meanKappa + sdKappa, 1)
                                    if Kappa < 0:
                                        print(Kappa)
                                # print("\nKappa: "+ str(Kappa) + " relHum: " + str(relHum) + " d: " + str(d) + " T: " + str(T))
                                d_wet = wetdiameter(kappa=Kappa, rh=relHum, ddry=d, temperature=T)
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
                                    fVoldry = wetVolumefraction(kappa=Kappa, rh=relHum, ddry=d, temperature=T)
                                    refIndAerosol = fVolBC * refIndBC + (1. - fVolBC) * refIndSol
                                    refIndWet: Union[complex, Any] = refIndWater * (
                                                1. - fVoldry) + refIndAerosol * fVoldry
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
                                result = ps.MieQ(refInd, diameter=dParticle, wavelength=wavelengths[countW],
                                                 asDict=True)

                            QBackVector[diameterBin, countW] = result['Qback'] * np.pi * 1.0 / 4.0 * (
                                    dParticle * 1e-9) ** 2.0 * numberConcentration / (4.0 * np.pi)
                            QAbsVector[diameterBin, countW] = result['Qabs'] * np.pi * 1.0 / 4.0 * (
                                    dParticle * 1e-9) ** 2.0 * numberConcentration
                            QScaVector[diameterBin, countW] = result['Qsca'] * np.pi * 1.0 / 4.0 * (
                                    dParticle * 1e-9) ** 2.0 * numberConcentration
                            QExtVector[diameterBin, countW] = result['Qext'] * np.pi * 1.0 / 4.0 * (
                                    dParticle * 1e-9) ** 2.0 * numberConcentration
                            # print(QExtVector)
                            drawProgressBar(n / numberOfCalculations, barLen=20)
                            n += 1
                            diameterBin += 1
                            # countW += 1

                        for nWavelengths in range(wavelengths.size):
                            SigmaAbs = integrateTrap(y=QAbsVector[:, nWavelengths], x=np.ma.transpose(diameter))
                            SigmaExt = integrateTrap(y=QExtVector[:, nWavelengths], x=np.ma.transpose(diameter))
                            SigmaSca = integrateTrap(y=QScaVector[:, nWavelengths], x=np.ma.transpose(diameter))
                            SigmaBack = integrateTrap(y=QBackVector[:, nWavelengths], x=np.ma.transpose(diameter))
                            SigmaBackVectorMonte[monte, nWavelengths] = SigmaBack
                            SigmaExtVectorMonte[monte, nWavelengths] = SigmaExt
                            # print("SigmaExtMonteMatrix Monte durchlauf" + str(monte))
                            # print(SigmaExtVectorMonte)
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
                # print(coefficients)
                finalCoefficients[i, 0] = timeOfScan
                finalCoefficientsSd[i, 0] = timeOfScan
                finalCoefficients[i, 1:] = coefficients
                finalCoefficientsSd[i, 1:] = coefficientsSd

            np.savetxt(fname=fout + "_SD.txt", X=finalCoefficientsSd, delimiter="\t")
            np.savetxt(fname=fout + ".txt", X=finalCoefficients, delimiter="\t")
            print("\nDuration: " + str(time.time() - start))

plt.plot(ScanTimes, finalCoefficients[:, 1], "o")
plt.plot(actosData[:, colTime], actosData[:, colRH], "bo")
plt.show()
