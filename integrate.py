import math
import numpy as np


def integrateTrap(x, y):
    #x = np.asarray(x)
    #y = np.asarray(y)
    print("y=", y)
    #print(y)
    print("x=", x)
    total = 0.0
    print(np.ma.shape(x)[0])
    for i in range(0, (np.ma.shape(x)[0] - 1), 1):
        print(i)
        total = total + (y[(i + 1)] + y[i]) / 2. * math.log10(x[(i + 1)] / x[i])
        return(total)
