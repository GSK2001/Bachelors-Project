from all_imports import *
#Opning HDU
hdu = fits.open("Q0918_f160w_NU.fits")[0]
noise=hdu.data[1000-75:1000+75,1000-50:1000+100].astype(float)

image=hdu.data[988-50:988+50,1055-70:1055+30].astype(float)
header=hdu.header
#Image scaling
z=ZScaleInterval()
z1,z2=z.get_limits(image)
#image origin
plt.rcParams['image.origin'] = 'lower'

#show image

plt.figure(1)
plt.imshow((image),origin="lower",cmap="binary",norm=LogNorm())
plt.colorbar()
plt.show()

# Estimate data noise at dark area
noise_cutout_pos = (128, 28)
noise_cutout_size = 60
noise_cutout = Cutout2D((noise), noise_cutout_pos, noise_cutout_size)

noise_mean = (noise_cutout.data).mean()
noise_sigma = (noise_cutout.data).std()
noise_3_sigma = noise_sigma * 3.
noise_8_sigma = noise_sigma * 8.

print(noise_mean, noise_3_sigma, noise_8_sigma)
# Plot image and noise distribution

plt.figure(2)
plt.imshow((noise_cutout.data),origin="lower",cmap="binary")
plt.title("Dark Patch")
plt.xlabel("Pixels")
plt.ylabel("Pixels")
plt.colorbar()
plt.show()

n, bins, patches = plt.hist(noise_cutout.data.flatten(), bins=35, align='left', color='black')
plt.plot(bins[:-1], n, c='r', linewidth=3)
plt.axvline(noise_mean, label="noise_mean", linestyle="--")

plt.xlabel('Flux Bins [log10(electrons/s)]')
plt.ylabel('Count')
plt.title('Noise Histogram')
plt.legend()
plt.show()

# Define detect threshold
threshold = noise_3_sigma
# Define smoothing kernel
kernel_size =3
fwhm = 5
npixels = 40
#Catalogue
cat, segm, segm_deblend = make_catalog(
    (image),
    threshold,
    deblend=True,
    kernel_size=kernel_size,
    fwhm=fwhm,
    npixels=npixels,vmin=z1, vmax=z2,
    plot=True
)
plt.show()

# Display source properties
print("Num of Targets:", len(cat))

# Convert to table
cat_table = cat.to_table()

print(cat_table[:10]) #to display the objects that is catalogued in a table

# Sort and get the largest object in the catalog
sorted_idx_list = order_cat(cat, key='area', reverse=True)
idx = sorted_idx_list[2] # index 0 is largest
source = cat[idx]  # get source from the catalog
#make a list of radius
r_list = make_radius_list(
    max_pix=20 # Max pixel to go up to
    ,n=200 # the number of radii to produce
)
#print(repr(r_list))
# Photomerty
flux_arr, area_arr, error_arr = source_photometry(

    # Inputs
    source, # Source (`photutils.segmentation.catalog.SourceCatalog`)
    image, # Image as 2D array
    segm_deblend, # Deblended segmentation map of image
    r_list, # list of aperture radii

    # Options
    cutout_size=max(r_list)*2, # Cutout out size, set to double the max radius
    bkg_sub=True, # Subtract background
    sigma=3, sigma_type='clip', # Fit a 2D plane to pixels within 1 sigma of the mean
    plot=True, vmax=z2, vmin=z1, # Show plot with max and min defined above
)
plt.show()

'''def max(x,a):
    return a

value, error = curve_fit(max,r_list[-20:],flux_arr[-20:])

plt.plot(r_list, flux_arr)
plt.errorbar(r_list,np.zeros(len(r_list))+max(r_list,value[0]),yerr=np.zeros(len(r_list))+error[0])
plt.errorbar(r_list,np.zeros(len(r_list))+max(r_list,value[0]/2),yerr=np.zeros(len(r_list))+error[0])
#
#interpolating the half light radius
#value=cat_table[2][18]
r_e=np.interp(value/2,flux_arr,r_list)

'''
plt.show()
#saves the photometric data as a table
t = Table(
    data=[r_list, flux_arr, area_arr, error_arr],
    names=['r_list', 'flux_arr', 'area_arr', 'error_arr'],
)

t.write('photometric data.csv', overwrite=True)

#print(t)
#construct a petrosian
p = Petrosian(r_list, area_arr, flux_arr)
print("half light in pixel is",p.r_half_light_arcsec(w(header))/0.06) 
conversion_pix_to_arcsec= 0.06 *u.arcsec
print("half light in arcsec is",p.r_half_light_arcsec(w(header)))#conversion_pix_to_arcsec*r_e )
z=2.583
V=Planck18.kpc_proper_per_arcmin(z)
ar=p.r_half_light_arcsec(w(header))*u.arcsec
print("half light in kpc is",V*(ar.to(u.arcmin)))
electron_flux=cat_table[2][16].tolist()
print("the b/a ratio is" ,cat_table[2][10]/cat_table[2][9])
c=flux_to_abmag( electron_flux, header ) #need to print this in shell for it to work
print("The electron flux from the DLA is", electron_flux , "and the area in pixel squared is", cat_table[2][7] )
print("The AB magnitude is", c)
