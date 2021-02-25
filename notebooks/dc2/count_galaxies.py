#!/usr/bin/env python
# coding: utf-8

# Modified from:
# Rubin LSST DESC DC2: Accessing Object Table with GCRCatalogs

###Import necessary packages
import numpy as np
import matplotlib.pyplot as plt
plt.ioff() 
#get_ipython().run_line_magic('matplotlib', 'inline')

###GCR catalogs import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
from GCRCatalogs import GCRQuery

###Access object table with GCRCatalogs
GCRCatalogs.get_public_catalog_names()
obj_cat = GCRCatalogs.load_catalog("dc2_object_run2.2i_dr6_v2")

###Object Table schema
sorted(obj_cat.list_all_quantities())

###Truth
extragalactic_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small')
#extragalactic_cat = GCRCatalogs.load_catalog('skysim5000_v1.1.1_small')
sorted(obj_cat.list_all_quantities())

###Define truth cuts
cosmoDC2_cuts = [
    #is_galaxy,
    #bright,
    GCRQuery((np.isfinite, 'mag_i')),  # Select objects that have i-band cmodel magnitudes
]

###Truth selection
columns = ["ra", "dec", "mag_i"]
truth = extragalactic_cat.get_quantities(
    quantities=columns,
    filters=cosmoDC2_cuts
)

###Data selection
is_extended = GCRQuery('extendedness == 1')  # Extended objects (primarily galaxies)
clean = GCRQuery('clean')  # The source has no flagged pixels (interpolated, saturated, edge, clipped...) 
                           # and was not skipped by the deblender

galaxy_cuts = [
    clean,
    is_extended, 
    GCRQuery((np.isfinite, 'mag_i')),  # Select objects that have i-band cmodel magnitudes
]

###Quantities to store
columns = ["ra", "dec", "mag_i", ]

data = obj_cat.get_quantities(
    quantities=columns, 
    #filters=galaxy_cuts, 
    native_filters=tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
)

data_cut = obj_cat.get_quantities(
    quantities=columns, 
    filters=galaxy_cuts, 
    native_filters=tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
)

###Compute sky area
d_ra = data["ra"].max() - data["ra"].min()
d_dec = data["dec"].max() - data["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(data["dec"])))
sky_area_sq_deg = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg)
print('********Plot********')
###Data plot
mag_bins = np.linspace(14, 24, 200)
cdf = np.searchsorted(data["mag_i"], mag_bins, sorter=data["mag_i"].argsort())
cdf_cut = np.searchsorted(data_cut["mag_i"], mag_bins, sorter=data_cut["mag_i"].argsort())
g1, = plt.semilogy(mag_bins, cdf / sky_area_sq_deg, color = 'black')
g2, = plt.semilogy(mag_bins, cdf_cut / sky_area_sq_deg, color = 'red')
plt.xlabel("$i$-band mag");
plt.ylabel("Galaxies /deg2 /0.05mag");
plt.legend([g1,g2],["dr6", "dr6 (ext=1)"], title = 'DC2 (9 tracts)')
#plt.grid();
#plt.grid(True, which="both", ls="-", )
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
###Define the output path for figures
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/"
plt.savefig(outpath+"dr6.png", bbox_inches='tight') 
plt.close()

###CosmoDC2 plot
# Cdf = Np.Searchsorted(Truth["Mag_i"], Mag_bins, Sorter=Truth["Mag_i"].Argsort())
# Sky_area_sq_deg = 50;
# Plt.Semilogy(Mag_bins, Cdf / Sky_area_sq_deg)
# Plt.Xlabel("$I$-Band Mag");
# Plt.Ylabel("Galaxies /Deg2 /0.05mag");
# Plt.Grid();  # Add Grid To Guide Our Eyes
# ###Define The Output Path For Figures
# Outpath = "/Sps/Lsst/Users/Tguillem/Desc/Config/Tests_dc2/Desc-Data-Portal/Notebooks/Dc2/Plots/"
# Plt.Savefig(Outpath+"Truth.Png", Bbox_inches='Tight')


###Filtering out the stars using the matched catalog
  ###to use the matched catalog to filter out the stars


###Redshift distributions
