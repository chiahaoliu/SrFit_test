#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

# We'll need numpy and matplotlib for plotting our results
import numpy as np
import matplotlib.pyplot as plt

# A least squares fitting algorithm from scipy
from scipy.optimize.minpack import leastsq, least_squares

# DiffPy-CMI modules for building a fitting recipe
from diffpy.structure import loadStructure
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe, FitResults

def makeRecipe_npload(x, y, struc, dy=None):
    # Files containing our experimental data and structure file
    #dataFile = "Ni_new.gr"
    #structureFile = "Ni_new.cif"
    spaceGroup = "Fm-3m"

    # The first thing to construct is a contribution. Since this is a simple
    # example, the contribution will simply contain our PDF data and an associated
    # structure file. We'll give it the name "nickel"
    niPDF = PDFContribution("nickel")
    niPDF.profile.setObservedProfile(x, y, dy)

    # Load the data and set the r-range over which we'll fit
    niPDF.setCalculationRange(xmin=1.5, xmax=30, dx=0.01)

    # Add the structure from our cif file to the contribution
    #niStructure = loadStructure(structureFile)
    niStructure = struc 
    niPDF.addStructure("nickel", niStructure)

    # The FitRecipe does the work of calculating the PDF with the fit variable
    # that we give it.
    niFit = FitRecipe()

    # give the PDFContribution to the FitRecipe
    niFit.addContribution(niPDF)

    # Configure the fit variables and give them to the recipe.  We can use the
    # srfit function constrainAsSpaceGroup to constrain the lattice and ADP
    # parameters according to the Fm-3m space group.
    from diffpy.srfit.structure import constrainAsSpaceGroup
    spaceGroupParams = constrainAsSpaceGroup(niPDF.nickel.phase, spaceGroup)
    print("Space group parameters are:",
          ', '.join(p.name for p in spaceGroupParams))
    print()

    # We can now cycle through the parameters and activate them in the recipe as
    
    print("  variables:", niFit.names)
    print("  initial values:", niFit.values)
    # variables
    for par in spaceGroupParams.latpars:
        niFit.addVar(par)
    # Set initial value for the ADP parameters, because CIF had no ADP data.
    for par in spaceGroupParams.adppars:
        niFit.addVar(par, value=0.005)

    # As usual, we add variables for the overall scale of the PDF and a delta2
    # parameter for correlated motion of neighboring atoms.
    niFit.addVar(niPDF.scale, 1)
    niFit.addVar(niPDF.nickel.delta1, 1)

    # We fix Qdamp based on prior information about our beamline.
    niFit.addVar(niPDF.qdamp, 0.04, fixed=False)
    niFit.addVar(niPDF.qbroad, 0.01, fixed=False)

    # Turn off printout of iteration number.
    niFit.clearFitHooks()

    return niFit

def makeRecipe_loaddata(obs_fn, struc):
    # Files containing our experimental data and structure file
    #dataFile = "Ni_new.gr"
    #structureFile = "Ni_new.cif"
    spaceGroup = "Fm-3m"

    # The first thing to construct is a contribution. Since this is a simple
    # example, the contribution will simply contain our PDF data and an associated
    # structure file. We'll give it the name "nickel"
    niPDF = PDFContribution("nickel")
    niPDF.loadData(obs_fn)

    # Load the data and set the r-range over which we'll fit
    niPDF.setCalculationRange(xmin=1.5, xmax=30, dx=0.01)

    # Add the structure from our cif file to the contribution
    #niStructure = loadStructure(structureFile)
    niStructure = struc 
    niPDF.addStructure("nickel", niStructure)

    # The FitRecipe does the work of calculating the PDF with the fit variable
    # that we give it.
    niFit = FitRecipe()

    # give the PDFContribution to the FitRecipe
    niFit.addContribution(niPDF)

    # Configure the fit variables and give them to the recipe.  We can use the
    # srfit function constrainAsSpaceGroup to constrain the lattice and ADP
    # parameters according to the Fm-3m space group.
    from diffpy.srfit.structure import constrainAsSpaceGroup
    spaceGroupParams = constrainAsSpaceGroup(niPDF.nickel.phase, spaceGroup)
    print("Space group parameters are:",
          ', '.join(p.name for p in spaceGroupParams))
    print()

    # We can now cycle through the parameters and activate them in the recipe as
    
    print("  variables:", niFit.names)
    print("  initial values:", niFit.values)
    # variables
    for par in spaceGroupParams.latpars:
        niFit.addVar(par)
    # Set initial value for the ADP parameters, because CIF had no ADP data.
    for par in spaceGroupParams.adppars:
        niFit.addVar(par, value=0.005)

    # As usual, we add variables for the overall scale of the PDF and a delta2
    # parameter for correlated motion of neighboring atoms.
    niFit.addVar(niPDF.scale, 1)
    niFit.addVar(niPDF.nickel.delta1, 1)

    # We fix Qdamp based on prior information about our beamline.
    niFit.addVar(niPDF.qdamp, 0.04, fixed=False)
    niFit.addVar(niPDF.qbroad, 0.01, fixed=False)

    # Turn off printout of iteration number.
    niFit.clearFitHooks()

    return niFit

fit_bound = (np.asarray([3.5, 0.003, 0., 0., 0.01, 0.005]),
             np.asarray([3.6, 0.01, 1., 4., 0.1, 0.1])
            )
"""
#leastsq(niFit.residual, niFit.values)
least_squares(niFit.residual, niFit.values, bounds=fit_bound)
print("  final values:", niFit.values)
print()

# Obtain and display the fit results.
niResults = FitResults(niFit)
print("FIT RESULTS\n")
print(niResults)

# Plot the observed and refined PDF.

# Get the experimental data from the recipe
r = niFit.nickel.profile.x
gobs = niFit.nickel.profile.y

# Get the calculated PDF and compute the difference between the calculated and
# measured PDF
gcalc = niFit.nickel.evaluate()
baseline = 1.1 * gobs.min()
gdiff = gobs - gcalc

# Plot!
plt.figure()
plt.plot(r, gobs, 'bo', label="G(r) data",
         markerfacecolor='none', markeredgecolor='b')
plt.plot(r, gcalc, 'r-', label="G(r) fit")
plt.plot(r, gdiff + baseline, 'g-', label="G(r) diff")
plt.plot(r, np.zeros_like(r) + baseline, 'k:')
plt.xlabel(r"r ($\AA$)")
plt.ylabel(r"G ($\AA^{-2}$)")
plt.legend()

plt.show()
"""
