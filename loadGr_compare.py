import numpy as np
from diffpy.structure import loadStructure
from diffpy.srfit.fitbase import FitResults
from scipy.optimize.minpack import leastsq, least_squares
from fitNi import makeRecipe_npload, makeRecipe_loaddata

r, Gr = np.loadtxt('Ni_new.gr', skiprows=25, dtype=np.float64).T
struc = loadStructure('Ni_new.cif')

# import makeRecipe with different loadData method
fit = makeRecipe_npload(r, Gr, struc)
fit2 = makeRecipe_loaddata('Ni_new.gr', struc)

leastsq(fit.residual, fit.values)
leastsq(fit2.residual, fit2.values)

for method, el in zip(['np.loadtxt', 'diffpy.srfit.pdf.pdfcontribution.loadData'],
                      [fit, fit2]):
    print("INFO: with method %s" % method)
    print("final param val = {}".format(el.values))
    print("rw = %.6f" % FitResults(el).rw)