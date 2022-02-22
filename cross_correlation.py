from tty import CC
from turtle import color
import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from scipy import interpolate

import pyratbay.constants as pc

try:
    wave , flux = np.loadtxt('Part1-reproducing_the_paper/piratbay_test/Giacobbe-main/HCN_Spectrum.dat', unpack=True)
except:
    wave , flux = np.loadtxt('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe/HCN_Spectrum.dat', unpack=True)

try:
    sys.path.append(os.path.abspath('/home/rahul/Desktop/notes/books_and_notes/liton_project/Part1-reproducing_the_paper/piratbay_test/Giacobbe-main'))
except:
    sys.path.append(os.path.abspath('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe'))

import wavelength_Grid_modelling as gridding

try:
    allspectrum1 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/allspectrum1.npy')
    wavelength1 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/wavelength1.npy')
    allspectrum2 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/allspectrum2.npy')
    wavelength2 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/wavelength2.npy')
    allspectrum3 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/allspectrum3.npy')
    wavelength3 = np.load('/home/rahul/Desktop/notes/books_and_notes/liton_project/wavelength3.npy')

except:
    allspectrum1 = np.load('/home/ws3/Desktop/rahul/allspectrum1.npy')
    wavelength1 = np.load('/home/ws3/Desktop/rahul/wavelength1.npy')
    allspectrum2 = np.load('/home/ws3/Desktop/rahul/allspectrum2.npy')
    wavelength2 = np.load('/home/ws3/Desktop/rahul/wavelength2.npy')
    allspectrum3 = np.load('/home/ws3/Desktop/rahul/allspectrum3.npy')
    wavelength3 = np.load('/home/ws3/Desktop/rahul/wavelength3.npy')
test_spectrum =[]


# rearranging
wav1=[]
flux1 =[]
for i in range(len(wave)-1,-1,-1):
    wav1.append(wave[i])
    flux1.append((flux[i]/pc.percent))

exclude = [8,9,10,23,24,30,45,46,47,48,49]

flux1 = np.array(flux1)
f = interpolate.splrep(wav1,flux1)



def CCF():
    CCmatrixf = np.zeros((len(np.transpose(allspectrum3[20])),167))
    for k in range(len(allspectrum3)):


        grid = gridding.grid_maker(wavelength3[k])
        CCmatrix = np.zeros((len(np.transpose(allspectrum3[k])),len(grid)))
        for i in range(len(CCmatrix)):
            for j in range(len(CCmatrix[0])):
                CCmatrix[i][j] = np.corrcoef(np.transpose(allspectrum3[k])[i], interpolate.splev(grid[j],f))[0][1]
        print(k)        

        CCmatrixf += CCmatrix

    plt.imshow(CCmatrixf, cmap='gray')
    plt.colorbar()
    plt.savefig('CCF100test.png')
    plt.clf()
    np.save('CCF100test.npy',CCmatrixf)

CCmatrixf = np.load('/home/ws3/Desktop/rahul/CCF.npy')

radial_vel = np.arange(100,-100,-1.5)
radial_vel2 = np.arange(225,-225,-2.7)
    



def realign_detection(CCmatrixf):

    kp = np.arange(0,200,3)
    phi = np.linspace(-0.03,0.03,len(radial_vel))

    detection = np.zeros((len(kp),len(radial_vel)))
    for i in range(len(kp)):
        rvpl = -14.741 + kp[i]*np.sin(2*np.pi*(0.007+phi))
        for j in range(len(CCmatrixf)):
            outrv = rvpl + radial_vel
            fit  = interpolate.interp1d(radial_vel2,CCmatrixf[j])
            detection[i] += fit(outrv)

    return detection    
    
CCF()

detection  = realign_detection(CCmatrixf)

plt.imshow(detection)
plt.colorbar()
plt.show()
    





