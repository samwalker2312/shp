import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.distributions import chi2
import math

params = {'figure.figsize': [8.3, 8], 'xtick.top':True, 'xtick.bottom':True,'ytick.left':True,'ytick.right':True}
mpl.rcParams.update(params)
path = 'output_cats/'
filename = 'mass_select.fits'
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

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True)
norm = mpl.colors.Normalize(vmin=-12.,vmax=-8)
mask1 = (z > 2) & (z < 2.5)
mask2 = (z > 2.5) & (z < 3)
mask3 = (z > 3) & (z < 3.5)
mask4 = (z > 3.5)

scatter1 = ax1.scatter(vminusj[mask1], uminusv[mask1], c=sSFR[mask1], cmap = 'RdYlGn', norm=norm, s = 9)
ax2.scatter(vminusj[mask2], uminusv[mask2],c=sSFR[mask2], cmap = 'RdYlGn', norm=norm, s = 9)
ax3.scatter(vminusj[mask3], uminusv[mask3],c=sSFR[mask3], cmap = 'RdYlGn', norm=norm, s = 9)
ax4.scatter(vminusj[mask4], uminusv[mask4],c=sSFR[mask4], cmap = 'RdYlGn', norm=norm, s = 9)

ax1.set_xlim(-.5,2.7)
ax1.set_ylim(-.5,2.5)
ax2.set_xlim(-.5,2.7)
ax2.set_ylim(-.5,2.5)
ax3.set_xlim(-.5,2.7)
ax3.set_ylim(-.5,2.5)
ax4.set_xlim(-.5,2.7)
ax4.set_ylim(-.5,2.5)

ax1.xaxis.set_major_locator(plt.MultipleLocator(.5))
ax2.xaxis.set_major_locator(plt.MultipleLocator(.5))
ax3.xaxis.set_major_locator(plt.MultipleLocator(.5))
ax4.xaxis.set_major_locator(plt.MultipleLocator(.5))

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
ax4.plot((-1.,xstart),(1.3,1.3), color="black")
ax4.plot(x,y, color="black")
ax4.plot((1.6,1.6),(1.6*.88 + zpt, 5.), color="black")

ax4.set_xlabel(r"$V - J$")
ax3.set_xlabel(r"$V - J$")

ax1.set_ylabel(r"$U - V$")
ax3.set_ylabel(r"$U - V$")

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()
ax4.label_outer()

size = 'small'
x1 = r'$\textrm{N} = ' + str(len(uminusv[mask1])) + r'$'
x2 = r'$\textrm{N} = ' + str(len(uminusv[mask2])) + r'$'
x3 = r'$\textrm{N} = ' + str(len(uminusv[mask3])) + r'$'
x4 = r'$\textrm{N} = ' + str(len(uminusv[mask4])) + r'$'
ax1.annotate(x1, xy = (-0.3,2.25), fontsize = size)
ax2.annotate(x2, xy = (-0.3,2.25), fontsize = size)
ax3.annotate(x3, xy = (-0.3,2.25), fontsize = size)
ax4.annotate(x4, xy = (-0.3,2.25), fontsize = size)

y1 = r'$2 < z < 2.5$'
y2 = r'$2.5 < z < 3$'
y3 = r'$3 < z < 3.5$'
y4 = r'$3.5 < z < 4.5$'

addsamplenumb = True
if addsamplenumb == True:
    with fits.open('output_cats/sample.fits') as hdu:
        sample_data = hdu[1].data
    z = sample_data['redshift_50']
    uv = sample_data['UV_colour_50']
    vj = sample_data['VJ_colour_50']
    mask1 = (z > 2) & (z < 2.5)
    mask2 = (z > 2.5) & (z < 3)
    mask3 = (z > 3) & (z < 3.5)
    mask4 = (z > 3.5)
    uvjmask1 = (uv[mask1] > 1.3) & (vj[mask1] < 1.6) & (uv[mask1] > (0.88*vj[mask1] + .69))
    uvjmask2 = (uv[mask2] > 1.3) & (vj[mask2] < 1.6) & (uv[mask2] > (0.88*vj[mask2] + .69))
    uvjmask3 = (uv[mask3] > 1.3) & (vj[mask3] < 1.6) & (uv[mask3] > (0.88*vj[mask3] + .69))
    uvjmask4 = (uv[mask4] > 1.3) & (vj[mask4] < 1.6) & (uv[mask4] > (0.88*vj[mask4] + .69))

    z1 = r'$\textrm{N}_{\textrm{sample}} = '+str(len(z[mask1])) + r'$' + '\n' + r'$\textrm{N}_{UVJ} = ' + str(len(z[mask1][uvjmask1])) + r'$'
    z2 = r'$\textrm{N}_{\textrm{sample}} = '+str(len(z[mask2])) + r'$' + '\n' + r'$\textrm{N}_{UVJ} = ' + str(len(z[mask2][uvjmask2])) + r'$'
    z3 = r'$\textrm{N}_{\textrm{sample}} = '+str(len(z[mask3])) + r'$' + '\n' + r'$\textrm{N}_{UVJ} = ' + str(len(z[mask3][uvjmask3])) + r'$'
    z4 = r'$\textrm{N}_{\textrm{sample}} = '+str(len(z[mask4])) + r'$' + '\n' + r'$\textrm{N}_{UVJ} = ' + str(len(z[mask4][uvjmask4])) + r'$'

    ax1.annotate(z1, xy = (-0.3,1.9), fontsize = size)
    ax2.annotate(z2, xy = (-0.3,1.9), fontsize = size)
    ax3.annotate(z3, xy = (-0.3,1.9), fontsize = size)
    ax4.annotate(z4, xy = (-0.3,1.9), fontsize = size)

    ax1.annotate(y1, xy = (-0.3,1.75), fontsize = size)
    ax2.annotate(y2, xy = (-0.3,1.75), fontsize = size)
    ax3.annotate(y3, xy = (-0.3,1.75), fontsize = size)
    ax4.annotate(y4, xy = (-0.3,1.75), fontsize = size)
else:
    ax1.annotate(y1, xy = (-0.3,2.), fontsize = size)
    ax2.annotate(y2, xy = (-0.3,2.), fontsize = size)
    ax3.annotate(y3, xy = (-0.3,2.), fontsize = size)
    ax4.annotate(y4, xy = (-0.3,2.), fontsize = size)

colorax = fig.add_axes([0.77, 0.13,0.03,0.17])
cbar = plt.colorbar(scatter1, cax=colorax)
cbar.ax.set_ylabel(r"$\mathrm{sSFR}$")
cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(1))

plt.subplots_adjust(wspace=0.08, hspace=0.08)
if apply_chisq == True:
    plt.savefig('plots/uvjplots_'+name+'_applychisq.pdf')
else:
    plt.savefig('plots/uvjplots_'+name+'.pdf')
#plt.show()
