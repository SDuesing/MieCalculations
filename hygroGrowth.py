# This module solves the semi-empirical single parameter representation of the hygroscopic growth and cloud
# condensation nuclei part with A from equation 8.69251x10-6 paper: https://bit.ly/2YuVZHx
from scipy.optimize import fsolve
import math
import numpy as np


def wetdiameter(ddry, kappa, temperature, rh):
    # ddry in µm
    # kappa dimensionless
    # temperature in °C
    # rh dimensionless (saturation 10%=0.1)
    def f(d):
        surfacetension = 0.072
        temp = temperature + 273.15
        return (d ** 3 - ddry ** 3) / (d ** 3 - (1 - kappa) * ddry ** 3) * math.exp(
            8.69251 * 10e-6 * surfacetension / d * 298.15 / temp) - rh
    x = fsolve(f, ddry + 0.1)  #ddry+0.1 is estimate starting value with 0.1 µm plus start diameter

    return x


def wetVolumefraction(ddry, rh, temperature, kappa):
    d_wet = wetdiameter(ddry, kappa, temperature, rh)
    f_vol = ddry ** 3. / d_wet ** 3.  # volume of dry particle to the total particle
    return f_vol
