from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

with fits.open('UDS-DR11-K.mef.fits') as hdu:
    hdu.info()
    #print(hdu[0].header)
    data = hdu[0].data

targetra = []
targetdec = []
IDs = [115779,140198,33726,10864,170965,164787,65302,200263,68768,58100,6651,107021]

dataset = pd.read_csv('K_selected_uds_photoz_masses_photom_upload.cat', header=0, delim_whitespace=True, index_col=0)

for ID in IDs:
    targetra.append(dataset.loc[ID,'RA'])
    targetdec.append(dataset.loc[ID,'Dec'])

#targetra = [34.38154]
#targetdec = [-5.09253]

pixelxplus = []
pixelyplus = []

pixelxminus = []
pixelyminus = []

w = wcs.WCS(hdu[0].header)
delta = 2.5/3600
for i in range(len(targetra)):
    pixcoordsplus = w.wcs_world2pix([[targetra[i]+delta,targetdec[i]+delta]],0)
    pixelxplus.append(int(pixcoordsplus[0,0]))
    pixelyplus.append(int(pixcoordsplus[0,1]))

    pixcoordsminus = w.wcs_world2pix([[targetra[i]-delta,targetdec[i]-delta]],0)
    pixelxminus.append(int(pixcoordsminus[0,0]))
    pixelyminus.append(int(pixcoordsminus[0,1]))

for i in range(len(pixelxplus)):
    pic = data[pixelxplus[i]:pixelxminus[i], pixelyminus[i]:pixelyplus[i]]
    #pic = data[100:200, 300:400]
    plt.figure()
    norm = mpl.colors.Normalize(vmin=0,vmax=50)
    plt.imshow(pic, cmap = 'gray', norm=norm)
    plt.colorbar()
    plt.show()
