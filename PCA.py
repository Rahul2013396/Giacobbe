from PIL.Image import new
from astropy.io import fits as ft
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from numpy.core.multiarray import nested_iters
plt.style.use(astropy_mpl_style)
import sqlite3

def datacollect(order,n):

    conn = sqlite3.connect('Part1-reproducing_the_paper/gofio-master/gofio3/webuidatabases/db_2019-09-03_offline.db')
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


def PCA(data,a):
    length = 2048
    obsspect = data
    # median normalisation
    for i in range(len(obsspect)):
        med = np.median(obsspect[i])
        obsspect[i]-= med

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

        for i in range(5): 
            S[i] = 0
        for i in range(len(S)):    
            s[i][i] = S[i]


        spect = np.dot(U,np.dot(s,VT))    
   


    newspect = np.transpose(spect-obsspect)
    vare =[]
    for i in range(len(newspect)):
        vare1 = np.var(newspect[i])
        newspect[i] /=vare1
        vare.append(vare1)

    newspect *= np.median(vare)
    newspect = np.transpose(newspect)
    
    return newspect


allspectrum = []
wavelength = [] 
for i in range(50):
    wavelength.append(datacollect(i,1)[0])
    data  = datacollect(i,2)
    spectrum = PCA(data,0)
    allspectrum.append(spectrum)
    print(i)

np.save('wavelength.npy', wavelength)
np.save('allspectrum.npy', allspectrum)    