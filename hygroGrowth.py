import math
from scipy.optimize import fsolve
import numpy as np


def wetdiameter(ddry, kappa, temperature, rh):
    def f(d):
        surfacetension = 0.072
        temp = temperature + 273.15
        return (d ** 3 - ddry ** 3) / (d ** 3 - (1 - kappa) * ddry ** 3) * math.exp(
            8.69251 * 10e-6 * surfacetension / d * 298.15 / temp) - rh

    x = fsolve(f, ddry+0.1)

    return x

