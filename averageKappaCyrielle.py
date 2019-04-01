import numpy as np
import datetime as dt
import pandas as pd


def averageKappa(table, datumsAngabe):

    datum = dt.datetime.strptime(datumsAngabe, "%Y%m%d")
    print("Datum :" + str(datum))
    tableKappaCyr = table
    #print(tableKappaCyr)
    kappa35 = tableKappaCyr[["date_d1_HTDMA", 'k35']].dropna()
    kappa50 = tableKappaCyr[["date_d50_HTDMA", 'k_HTDMA_D50']].dropna()
    kappa75 = tableKappaCyr[["date_d75_HTDMA", 'k_HTDMA_D75']].dropna()
    kappa110 = tableKappaCyr[["date_d110_HTDMA", 'k_HTDMA_D110']].dropna()
    kappa165 = tableKappaCyr[["date_d165_HTDMA", 'k_HTDMA_D165']].dropna()
    kappa265 = tableKappaCyr[["date_d265_HTDMA", 'k_HTDMA_D265']].dropna()

    #print(kappa35)

    datesKappa35 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa35["date_d1_HTDMA"]]).date

    #print(datesKappa35)

    datesKappa50 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa50["date_d50_HTDMA"]]).date
    datesKappa75 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa75["date_d75_HTDMA"]]).date
    datesKappa110 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa110["date_d110_HTDMA"]]).date
    datesKappa165 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa165["date_d165_HTDMA"]]).date
    datesKappa265 = pd.DatetimeIndex(
        [dt.datetime.strptime(date, "%d.%m.%Y %H:%M:%S") for date in
         kappa265["date_d265_HTDMA"]]).date

    k35 = kappa35["k35"]
    #print(datum)
    indize35Start = int(np.argwhere(datesKappa35 == datum.date())[0])
    indize35End = int(np.argwhere(datesKappa35 == datum.date())[-1])
    k50 = kappa50["k_HTDMA_D50"]
    indize50Start = int(np.argwhere(datesKappa50 == datum.date())[0])
    indize50End = int(np.argwhere(datesKappa50 == datum.date())[-1])
    k75 = kappa75["k_HTDMA_D75"]
    indize75Start = int(np.argwhere(datesKappa75 == datum.date())[0])
    indize75End = int(np.argwhere(datesKappa75 == datum.date())[-1])
    k110 = kappa110["k_HTDMA_D110"]
    indize110Start = int(np.argwhere(datesKappa110 == datum.date())[0])
    indize110End = int(np.argwhere(datesKappa110 == datum.date())[-1])
    k165 = kappa165["k_HTDMA_D165"]
    indize165Start = int(np.argwhere(datesKappa165 == datum.date())[0])
    indize165End = int(np.argwhere(datesKappa165 == datum.date())[-1])
    k265 = kappa265["k_HTDMA_D265"]
    indize265Start = int(np.argwhere(datesKappa265 == datum.date())[0])
    indize265End = int(np.argwhere(datesKappa265 == datum.date())[-1])


    #print(np.array(k35))

    return(np.array([np.mean(k35[indize35Start:(indize35End + 1)]), np.mean(k50[indize50Start:(indize50End + 1)]),
            np.mean(k75[indize75Start:(indize75End + 1)]), np.mean(k110[indize110Start:(indize110End + 1)]),
            np.mean(k165[indize165Start:(indize165End + 1)]), np.mean(k265[indize265Start:(indize265End + 1)]), np.std(k35[indize35Start:(indize35End + 1)]), np.std(k50[indize50Start:(indize50End + 1)]),
            np.std(k75[indize75Start:(indize75End + 1)]), np.std(k110[indize110Start:(indize110End + 1)]),
            np.std(k165[indize165Start:(indize165End + 1)]), np.std(k265[indize265Start:(indize265End + 1)])]))


def findKappa(d, diameter, kappas):
    dia = d
    if d < diameter[0]:
        kappa = kappas[0]
    elif d > diameter[-1]:
        kappa = kappas[-1]
    else:
        indexStart = int(np.argwhere(diameter < dia)[-1])
        indexEnd = int(np.argwhere(diameter > dia)[0])

        kappa = (kappas[indexEnd]-kappas[indexStart])/(diameter[indexEnd]-diameter[indexStart])*(dia-diameter[indexStart])+kappas[indexStart]

        #print(kappa)

    return kappa


