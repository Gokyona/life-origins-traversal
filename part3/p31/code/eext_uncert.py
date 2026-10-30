"""Uncertainty of the apparent extension activation energy from Aitken et al. 2023 Fig. 7A/B (pre-activated *N 10 mM,
*NNN 500 uM; >= +7 at 24 h: 42% at RT, 8% at +3 C). Homogeneous Poisson per-letter model: Lambda(T) from the yield,
E = -R ln(Lambda_3 / Lambda_RT) / (1/T3 - 1/T_RT). Uncertainty: yields +-5 points (gel quantification), RT 20-25 C.
Also an alternative read-out: the mean number of letters added from the gel band pattern is not available
numerically, so only the >= +7 fraction is used. Output eext_uncert.json."""
import json, math, itertools
import numpy as np
from scipy.stats import poisson
R = 8.314
lam = lambda p: float(next(l for l in np.linspace(0.001, 40, 400001) if poisson.sf(6, l) >= p))
def E(pRT, p3, TRT=22.0, T3=3.0):
    return -R * math.log(lam(p3) / lam(pRT)) / (1 / (T3 + 273.15) - 1 / (TRT + 273.15)) / 1000
nom = E(0.42, 0.08)
grid = [E(a, b, t) for a, b, t in itertools.product((0.37, 0.42, 0.47), (0.03, 0.08, 0.13), (20.0, 22.0, 25.0))]
out = {"E_nominal_kJ": nom, "E_range_yield_pm5_RT20to25_kJ": [min(grid), max(grid)],
       "E_yield_only_pm5_kJ": [E(0.47, 0.03), E(0.37, 0.13)]}
json.dump(out, open("eext_uncert.json", "w"), indent=1)
print(out)
