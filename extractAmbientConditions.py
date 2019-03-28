import numpy as np


def extractAmbient(table, rhCol, timeCol, tempCol, startTime, duration):
    time = table[:, int(timeCol)]
    indizestart = int(np.argwhere((time >= startTime) & (time < startTime + duration))[0])
    indizeend = int(np.argwhere((time >= startTime) & (time < startTime + duration))[-1])

    rh = table[indizestart:(indizeend + 1), rhCol]
    temp = table[indizestart:(indizeend + 1), tempCol]

    meanrh = np.mean(rh)
    sdrh = np.std(rh)
    meant = np.mean(temp)
    sdt = np.std(temp)

    return [meanrh, sdrh, meant, sdt]

