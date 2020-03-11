import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.distributions import chi2
import math

path = 'output_cats/'
filename = 'subsample.fits'
filepath = path+filename
name = filename[:-5]
with fits.open(filepath) as hdu:
    data = hdu[1].data

uminusv = data['UV_colour_50']
vminusj = data['VJ_colour_50']
z = data['redshift_50']
sSFR = data['sSFR_50']
chisq = data['chisq_phot']
nu = data['n_bands'] - 10

apply_chisq = True
if apply_chisq == True:
    mask = chi2.sf(chisq, nu) > 0.5*math.erfc(5*(2**-0.5))
    uminusv = uminusv[mask]
    vminusj = vminusj[mask]
    z = z[mask]
    sSFR = sSFR[mask]

fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True)
norm = mpl.colors.Normalize(vmin=-12.,vmax=-8)
mask1 = (z > 2) & (z < 3)
mask2 = (z > 3) & (z < 4)
mask3 = z > 4

scatter1 = ax1.scatter(vminusj[mask1], uminusv[mask1], c=sSFR[mask1], cmap = 'RdYlGn', norm=norm, s = 9)
ax2.scatter(vminusj[mask2], uminusv[mask2],c=sSFR[mask2], cmap = 'RdYlGn', norm=norm, s = 9)
ax3.scatter(vminusj[mask3], uminusv[mask3],c=sSFR[mask3], cmap = 'RdYlGn', norm=norm, s = 9)

ax1.set_xlim(-.5,2.5)
ax1.set_ylim(-.5,2.5)
ax2.set_xlim(-.5,2.5)
ax2.set_ylim(-.5,2.5)
ax3.set_xlim(-.5,2.5)
ax3.set_ylim(-.5,2.5)

zpt = .69
xstart = (1.3-zpt)/.88
xend = 1.6
ax1.plot((-1.,xstart),(1.3,1.3), color="black")
x = np.linspace(xstart, xend)
y = .88*x + zpt
ax1.plot(x,y, color="black")
ax1.plot((1.6,1.6),(1.6*.88 + zpt, 5.), color="black")
ax2.plot((-1.,xstart),(1.3,1.3), color="black")
ax2.plot(x,y, color="black")
ax2.plot((1.6,1.6),(1.6*.88 + zpt, 5.), color="black")
ax3.plot((-1.,xstart),(1.3,1.3), color="black")
ax3.plot(x,y, color="black")
ax3.plot((1.6,1.6),(1.6*.88 + zpt, 5.), color="black")

ax1.set_xlabel(r"$V - J$")
ax2.set_xlabel(r"$V - J$")
ax3.set_xlabel(r"$V - J$")
ax1.set_ylabel(r"$U - V$")

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()

size = 'medium'
x1 = r'$N = ' + str(len(uminusv[mask1])) + r'$'
x2 = r'$N = ' + str(len(uminusv[mask2])) + r'$'
x3 = r'$N = ' + str(len(uminusv[mask3])) + r'$'
ax1.annotate(x1, xy = (-0.3,2.25), fontsize = size)
ax2.annotate(x2, xy = (-0.3,2.25), fontsize = size)
ax3.annotate(x3, xy = (-0.3,2.25), fontsize = size)

y1 = r'$2 < z < 3$'
y2 = r'$3 < z < 4$'
y3 = r'$4 < z$'

ax1.annotate(y1, xy = (-0.3,2.), fontsize = size)
ax2.annotate(y2, xy = (-0.3,2.), fontsize = size)
ax3.annotate(y3, xy = (-0.3,2.), fontsize = size)

cbar = plt.colorbar(scatter1)
cbar.ax.set_ylabel(r"$\log_{10} \mathrm{(sSFR)}$")

plt.subplots_adjust(wspace=0, hspace=0)
if apply_chisq == True:
    plt.savefig('plots/uvjplots_'+name+'_applychisq.pdf')
else:
    plt.savefig('plots/uvjplots_'+name+'.pdf')
#plt.show()
