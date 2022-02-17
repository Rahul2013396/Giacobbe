from tty import CC
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from scipy import interpolate


wave , flux = np.loadtxt('Part1-reproducing_the_paper/piratbay_test/Giacobbe-main/HCN_Spectrum.dat', unpack=True)

sys.path.append(os.path.abspath('/home/rahul/Desktop/notes/books_and_notes/liton_project/Part1-reproducing_the_paper/piratbay_test/Giacobbe-main'))

import wavelength_Grid_modelling as gridding

allspectrum = np.load('Part1-reproducing_the_paper/piratbay_test/Giacobbe-main/allspectrum.npy')
wavelength = np.load('Part1-reproducing_the_paper/piratbay_test/Giacobbe-main/wavelength.npy')
wav1=[]
flux1 =[]
for i in range(len(wave)-1,-1,-1):
    wav1.append(wave[i])
    flux1.append(flux[i])
print(wav1)

f = interpolate.splrep(wav1,flux1)
CCmatrixf = np.zeros((len(np.transpose(allspectrum[29])),167))
for i in range(1):
    grid = gridding.grid_maker(wavelength[29]*1e-3)
    print(np.shape(allspectrum))
    CCmatrix = np.zeros((len(np.transpose(allspectrum[29])),len(grid)))

    for i in range(len(CCmatrix)):
        for j in range(len(CCmatrix[0])):
            CCmatrix[i][j] = np.correlate(np.transpose(allspectrum[29])[i], interpolate.splev(grid[j],f))
            #print(np.shape(allspectrum[42][i]),np.shape(f(grid[j])))
            #print(np.correlate(np.transpose(allspectrum[42])[i], f(grid[j])),'valid')

    for i in range(len(CCmatrix)):
        CCmatrix[i] /= np.sum(CCmatrix[i])
    CCmatrixf += CCmatrix

x = np.arange(1,168)
plt.plot(x,CCmatrixf[21])
plt.show()
plt.imshow(CCmatrixf, cmap='gray')
plt.show()