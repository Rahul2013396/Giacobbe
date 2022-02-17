from tty import CC
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d


wave , flux = np.loadtxt('pyrat_bay/Giacobbe_codes/HCN_Spectrum.dat', unpack=True)

sys.path.append(os.path.abspath('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe_codes'))

import wavelength_Grid_modelling as gridding

allspectrum = np.load('pyrat_bay/Giacobbe_codes/allspectrum.npy')
wavelength = np.load('pyrat_bay/Giacobbe_codes/wavelength.npy')



f = interp1d(wave,flux)
CCmatrixf = np.zeros((len(allspectrum[42]),167))
for i in range(1):
    grid = gridding.grid_maker(wavelength[42]*1e-3)

    CCmatrix = np.zeros((len(allspectrum[42]),len(grid)))

    for i in range(len(CCmatrix)):
        for j in range(len(CCmatrix[0])):
            CCmatrix[i][j] = np.correlate(allspectrum[42][i], f(grid[j]))
            print(np.correlate(allspectrum[42][i], f(grid[j])))

    #for i in range(len(CCmatrix)):
    #    CCmatrix[i] -= np.median(CCmatrix[i])
    CCmatrixf += CCmatrix

x = np.arange(1,168)
plt.plot(x,CCmatrixf[21])
plt.show()
plt.imshow(CCmatrixf, cmap='gray')
plt.show()