from tty import CC
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from scipy import interpolate

try:
    wave , flux = np.loadtxt('Part1-reproducing_the_paper/piratbay_test/Giacobbe-main/HCN_Spectrumemm.dat', unpack=True)
except:
    wave , flux = np.loadtxt('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe/HCN_Spectrum.dat', unpack=True)

try:
    sys.path.append(os.path.abspath('/home/rahul/Desktop/notes/books_and_notes/liton_project/Part1-reproducing_the_paper/piratbay_test/Giacobbe-main'))
except:
    sys.path.append(os.path.abspath('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe'))

import wavelength_Grid_modelling as gridding

try:
    allspectrum = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/allspectrum.npy')
    wavelength = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/wavelength.npy')

except:
    allspectrum = np.load('/home/ws3/Desktop/rahul/allspectrum.npy')
    wavelength = np.load('/home/ws3/Desktop/rahul/wavelength.npy')
test_spectrum =[]


wav1=[]
flux1 =[]
for i in range(len(wave)-1,-1,-1):
    wav1.append(wave[i])
    flux1.append(flux[i])


flux1 = np.array(flux1)
f = interpolate.splrep(wav1,flux1)

print(wavelength[0][0] ,wavelength[-1][-1])

plt.plot(wav1,flux1)
plt.plot(wavelength[40]*1e-3, np.transpose(allspectrum[40])[25])
plt.show()


CCmatrixf = np.zeros((len(np.transpose(allspectrum[20])),167))
for k in range(21):
    grid = gridding.grid_maker(wavelength[k]*1e-3)
    CCmatrix = np.zeros((len(np.transpose(allspectrum[k])),len(grid)))

    for i in range(len(CCmatrix)):
        for j in range(len(CCmatrix[0])):
            CCmatrix[i][j] = np.corrcoef(np.transpose(allspectrum[k])[i], interpolate.splev(grid[j],f))[0][1]
            

    for i in range(len(CCmatrix)):
        CCmatrix[i] -= np.median(CCmatrix[i])
    CCmatrixf += CCmatrix

 

plt.imshow(CCmatrixf[:,:], cmap='gray')
plt.colorbar()
plt.savefig('testful.png')
plt.show()


x = np.arange(1,168)
for i in range(len(CCmatrixf)):
    a = max(CCmatrixf[i])
    for j in range(len(CCmatrix[i])):
        if(CCmatrixf[i][j] != a):
            CCmatrixf[i][j] =0

plt.imshow(CCmatrixf[:,:], cmap='gray')
plt.colorbar()
plt.savefig('test.png')
plt.show()


test =[]
CCmatrixft = np.zeros((len(np.transpose(allspectrum[20])),167))
for i in range(1):
    grid = gridding.grid_maker(wavelength[20]*1e-3)
    for i in range(50):
        test_spectrum.append(interpolate.splev(grid[40+i],f))
    #plt.imshow(test_spectrum, cmap='gray')
    #plt.show()
    CCmatrixt = np.zeros((len(test_spectrum),len(grid)))

    for i in range(len(CCmatrixt)):
        for j in range(len(CCmatrixt[0])):
            CCmatrixt[i][j] = np.corrcoef(test_spectrum[i], interpolate.splev(grid[j],f))[0][1]
            

    CCmatrixft += CCmatrixt

#for i in range(len(CCmatrixft)):
#    a = max(CCmatrixft[i])
#    for j in range(len(CCmatrixt[i])):
#        if(CCmatrixft[i][j] != a):
#            CCmatrixft[i][j] =0
#        

#plt.imshow(CCmatrixft, cmap='gray')
#plt.show()