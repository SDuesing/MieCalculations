import time
from typing import Tuple

from numpy.core._multiarray_umath import ndarray
from progressbar import drawProgressBar
import matplotlib.pyplot as plt

from hygroGrowth import wetdiameter
import numpy as np
vfunc = np.vectorize(wetdiameter)

testN = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
testResult = np.zeros((testN.shape[0], 3))
for n in range(testN.size):
    d = np.random.randint(100, 301, size=testN[n])
    start = time.time()
    res = vfunc(d, kappa=0.4, rh=0.8, temperature=20)
    end = time.time()
    time1 = end - start

    start = time.time()
    d_wet = []
    for i in range(d.size):
        res = wetdiameter(ddry=d[i], kappa=0.4, rh=0.8, temperature=20)
        d_wet.append(res)
    end = time.time()
    time2 = end - start
    testResult[n, 0:3] = np.transpose([testN[n], time1, time2])
    drawProgressBar((n+1.)/testN.size, barLen=20)

plt.plot(testResult[:, 0], testResult[:, 1], "bo")
plt.plot(testResult[:, 0], testResult[:, 2], 'ro')
plt.show()

print("\n" + str(testResult[testResult.shape[0]-1, 1]/testResult[testResult.shape[0]-1, 2]))
