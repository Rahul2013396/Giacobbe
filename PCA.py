import cmath
from PIL.Image import new
from astropy.io import fits as ft
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from numpy.core.multiarray import nested_iters
plt.style.use(astropy_mpl_style)
import sqlite3
import sys
import os
from scipy import interpolate



sys.path.append(os.path.abspath('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe'))

import wavelength_Grid_modelling as gridding

wave , flux = np.loadtxt('/home/ws3/Desktop/rahul/pyrat_bay/Giacobbe/HCN_Spectrumlow.dat', unpack=True)
test =[]
wav1=[]
flux1 =[]
for i in range(len(wave)-1,-1,-1):
    wav1.append(wave[i])
    flux1.append((flux[i]/0.01))



f = interpolate.splrep(wav1,flux1)






def datacollect(order,n):
    try :
        conn = sqlite3.connect('/home/ws3/Desktop/rahul/gofio/gofio/gofio3/webuidatabases/db_2019-09-03_offline.db')
    except:
        conn = sqlite3.connect('gofio/gofio/gofio3/webuidatabases/db_2018-07-07_offline.db')

    cur = conn.cursor()
    cur.execute("SELECT * FROM spec1dfiles  WHERE spec_1d_type='ms1d' AND slitpos = 'A' ")
    
    rows = cur.fetchall()

    spectrumA = []
    spectrumB = []
    for i in range(len(rows)):
        a = ft.getdata(rows[i][0])
        spectrumA.append(a[order][n])


    cur = conn.cursor()
    cur.execute("SELECT * FROM spec1dfiles  WHERE spec_1d_type='ms1d' AND slitpos = 'B' ")

    
    rows = cur.fetchall()

    for i in range(len(rows)):
        a = ft.getdata(rows[i][0])
        spectrumB.append(a[order][n])

    spectrum = []
    for i in range(len(spectrumB)):
        spectrum.append(spectrumA[i])
        spectrum.append(spectrumB[i])     

    return spectrum

exclude = [8,9,10,23,24,30,45,46,47,48,49]


fig , axs = plt.subplots(5)

def PCA(data,a):
    length = 2048
    obsspect = data
    axs[0].imshow(obsspect[:,1250:1750],cmap = 'gray')
    plt.grid(False)
    
    # median normalisation
    for i in range(len(obsspect)):
        med = np.median(obsspect[i])
        obsspect[i]/= med

    axs[1].imshow(obsspect[:,1250:1750],cmap = 'gray')
    plt.grid(False)

    #standardization(mean)
    obsspect = np.transpose(obsspect)

    for i in range(len(obsspect)):
        mean = np.mean(obsspect[i])
        obsspect[i] -= mean

    obsspect = np.transpose(obsspect)


    #standardization(std)
    for i in range(len(obsspect)):
        std = np.std(obsspect[i])
        obsspect[i] /= std
    axs[2].imshow(obsspect[:,1250:1750],cmap = 'gray')
    plt.grid(False)

    if(a==1):
        F = np.ones((len(data),2048),dtype='float')
        covm = np.cov(np.transpose(obsspect))


        u, v  = np.linalg.eig(covm)
        v= np.real(v)
        u= np.real(u)
    
        for i in range(len(data)):
            F[i] = v[:,i] 

        obsspect =np.transpose(obsspect) 
        a = np.dot(F,obsspect)

        G = np.zeros((len(data),length),dtype='float')
        for i in range(5):
             G[i] = v[:,i]

        spect  = np.dot(np.transpose(G),a)
        obsspect = np.transpose(obsspect)
        spect = np.transpose(spect)

    else:
        obsspect = np.transpose(obsspect)
        U,S,VT = np.linalg.svd(obsspect)
    
        s=np.zeros((len(U),len(S)))

        for i in range(4): 
            S[i] = 0
        for i in range(len(S)):    
            s[i][i] = S[i]


        spect = np.dot(U,np.dot(s,VT))    
   
    axs[3].imshow(np.transpose(spect)[:,1250:1750],cmap = 'gray')
    plt.grid(False)

    axs[4].imshow(np.transpose(obsspect-spect)[:,1250:1750],cmap = 'gray')
    plt.grid(False)
    newspect = spect
    
    plt.savefig('PCA.png')

    vare =[]
    for i in range(len(newspect)):
        vare1 = np.var(newspect[i])
        newspect[i] /=vare1
        vare.append(vare1)

    newspect *= np.median(vare)
    return newspect


allspectrum = []
wavelength = [] 
for i in range(50):
    if (i not in exclude):
        wavelength.append(datacollect(i,1)[0]*1e-3)
        data  = datacollect(i,2)
        test1 =np.array(data)
        grid = gridding.grid_maker(wavelength[-1])
        test = []
        for j in range(0,2*len(data),2):
            test.append(interpolate.splev(grid[j],f)) 
        test = np.array(test)
        data += 0*test
        data = np.array(data)
        spectrum = PCA(data,0)
        
        allspectrum.append(spectrum)
        print(i)

np.save('wavelength3.npy', wavelength)
np.save('allspectrum3.npy', allspectrum)    