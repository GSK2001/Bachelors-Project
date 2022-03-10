from all_imports import *

def FITS(file):
    image=FITSFigure(file,hdu=0,north=True)
    #image.recenter(image.pixel2world(620, 480, wcs=None)[0],image.pixel2world(620, 480, wcs=None)[1],radius=0.002)
    image.show_grayscale()
    image.add_scalebar(0.1)
    image.add_grid()
    #image.show_contour()
    image.add_colorbar()
    image.save("Q0918_f160w_NU.jpg")
    return plt.show()

#FITS("Q0918_f160w_NU.fits")


#opening image for peak identification
hdu = fits.open("Q0918_f160w_NU.fits")[0]
image=hdu.data[1000-75:1000+75,1000-75:1000+75].astype(float)
header= w(hdu.header)




#peak identity
data=image
mean, median, std = scs(data, sigma=15.0)
threshold = median + (6. * std)
#print(mean, median, std)
tbl = find_peaks(data, threshold, box_size=11)
tbl['peak_value'].info.format = '%.8g'
positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=1.)
norm = simple_norm(data, 'log', percent=99.9)
apertures.plot(color='#0547f9', lw=1.5)
plt.imshow(data, cmap='binary', origin='lower', norm=LogNorm(vmin=0.001, vmax=np.max(data)))
cbar = plt.colorbar()
cbar.set_label('Log Brightness', rotation=270, labelpad=25)
print(tbl[:10])# table of sources in terms of the pixel coordinate
plt.show()


#separation
c1=SC(header.pixel_to_world(1000-75+tbl[0][0],1000-50+tbl[0][1]))
c2=SC(header.pixel_to_world(1000-75+tbl[1][0],1000-50+tbl[1][1]))
sep=c1.separation(c2)
#
z     = 2.587                             # redshift of galaxies
theta = sep                                 # angle
r_ang = Planck18.kpc_proper_per_arcmin(z) # phys. dist. per angle
r     = r_ang * theta                     # physical distance
print('impact perameter between the DLA and QSO: ', r.to(u.kpc),sep.arcsec, "arcsec at z=2.587")
