# This module solves the semi-empirical single parameter representation of the hygroscopic growth and cloud
# condensation nuclei part with A from equation 8.69251x10-6 paper: https://bit.ly/2YuVZHx
from scipy.optimize import fsolve
import math
import numpy as np


def wetdiameter(ddry, kappa, temperature, rh):
    # ddry in nm
    # kappa dimensionless
    # temperature in °C
    # rh dimensionless (saturation 10%=0.1)
    ddry = ddry
    def f(d):
        surfacetension = 0.072
        temp = temperature + 273.15
        return (d ** 3 - ddry ** 3) / (d ** 3 - (1 - kappa) * ddry ** 3) * math.exp(
            8.69251 * 10e-6 * surfacetension / d * 298.15 / temp) - rh

    x = fsolve(f, ddry)
    if kappa < 0:
        print("Dry: " + str(ddry) + "wet: " + str(x) + "\n")
    return x



def wetVolumefraction(ddry, rh, temperature, kappa):
    d_wet = wetdiameter(ddry, kappa, temperature, rh)
    f_vol = ddry ** 3. / d_wet ** 3.  # volume of dry particle to the total particle
    return f_vol

"""
count = 0
kappa = np.linspace(0.1, 0.4, 40)
diameter = np.linspace(0.1, 5.0, 500)
result = np.zeros((diameter.shape[0] * kappa.shape[0], 3))
for k in kappa:
    for d in diameter:
        d_wet = wetdiameter(d, k, 20, 0.9) / d
        row = np.array([k, d, d_wet])
        result[count, 0:3] = np.transpose(row)
        count += 1
np.savetxt(fname='d_wet_table.txt', X=result, delimiter="\t")
"""