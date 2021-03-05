#!/usr/bin/env python
# coding: utf-8

# Modified from:
# Rubin LSST DESC DC2: Accessing Object Table with GCRCatalogs

###Import necessary packages
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.ioff() 
#get_ipython().run_line_magic('matplotlib', 'inline')

###lsst
#import lsst.afw.geom as afw_geom

###GCR catalogs import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
from GCRCatalogs import GCRQuery

###Access object table with GCRCatalogs
GCRCatalogs.get_public_catalog_names()
print(GCRCatalogs.get_available_catalog_names(include_default_only=False))
obj_cat = GCRCatalogs.load_catalog("dc2_object_run2.2i_dr6_v2")

###Object Table schema
#print('obj_cat quantities')
#print(sorted(obj_cat.list_all_quantities()))

###Truth
extragalactic_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small')
#extragalactic_cat = GCRCatalogs.load_catalog('skysim5000_v1.1.1_small')
#extragalactic_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_summary_tract_partition')
print('extragalactic_cat quantities')
print(sorted(extragalactic_cat.list_all_quantities()))

###Truth catalog
#truth_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_star_truth_summary')
truth_cat = GCRCatalogs.load_catalog('dc2_truth_run2.2i_summary_tract_partition')
#print('truth_cat quantities')
#print(sorted(truth_cat.list_all_quantities()))
###Dump variable definitions
for qty in ['cosmodc2_hp', 'cosmodc2_id', 'dec', 'flux_g', 'flux_g_noMW', 'flux_i', 'flux_i_noMW', 'flux_r', 'flux_r_noMW', 'flux_u', 'flux_u_noMW', 'flux_y', 'flux_y_noMW', 'flux_z', 'flux_z_noMW', 'host_galaxy', 'id', 'is_pointsource', 'is_variable', 'mag_g', 'mag_g_noMW', 'mag_i', 'mag_i_noMW', 'mag_r', 'mag_r_noMW', 'mag_u', 'mag_u_noMW', 'mag_y', 'mag_y_noMW', 'mag_z', 'mag_z_noMW', 'patch', 'ra', 'redshift', 'tract', 'truth_type']:
    info_dict = truth_cat.get_quantity_info(qty)
    print(qty,info_dict)
    #print('loop') 

###Truth-match catalog
truth_match_cat = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr6_v2_with_addons')
#print('truth_match_cat quantities')
#print(sorted(truth_match_cat.list_all_quantities()))

###Photo-z catalog
photoz_cat = GCRCatalogs.load_catalog('dc2_object_run2.2i_dr3a_with_photoz')
print('photoz_cat quantities')
print(sorted(q for q in photoz_cat.list_all_quantities() if q.startswith('photoz_')))

###Define truth cuts
cosmoDC2_cuts = [
    #is_galaxy,
    #bright,
    GCRQuery('mag_i<24')  # Select objects that have i-band cmodel magnitudes
]

###Truth selection
columns_truth = ["ra", "dec", "mag_i", "redshift"]
#truth = extragalactic_cat.get_quantities(
#    quantities=columns_truth,
#    filters=cosmoDC2_cuts
#)

###Data selection
is_extended = GCRQuery('extendedness == 1')  # Extended objects (primarily galaxies)
clean = GCRQuery('clean')  # The source has no flagged pixels (interpolated, saturated, edge, clipped...) 
                           # and was not skipped by the deblender
true_galaxy = GCRQuery('truth_type == 1')
true_star = GCRQuery('truth_type == 2')
#true_SN = [GCRQuery('truth_type !=1'),GCRQuery('truth_type !=2')]
is_good_match = GCRQuery('is_good_match == 1')
true_galaxy_cosmoDC2 = GCRQuery('truth_id == 1')

galaxy_cuts = [
    clean,#is_good_match,
    is_extended, 
    #true_galaxy,
    #GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
    GCRQuery('mag_i_cModel<24')
]

star_cuts = [
    clean,is_good_match,
    #is_extended, 
    true_star,
    GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
]

is_extended_cuts = [
    clean,is_good_match,
    #is_extended, 
    #true_SN,
    GCRQuery('truth_type!=1'),
    #GCRQuery('truth_type !=2'),
    GCRQuery((np.isfinite, 'mag_i_cModel')),  # Select objects that have i-band cmodel magnitudes
]

###Quantities to store
columns = ["ra", "dec", "mag_i_cModel", "photoz_mode"]

#data = photoz_cat.get_quantities(
#    quantities=columns, 
#    native_filters=tract_filter([3830])#, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

data_galaxy_cut = photoz_cat.get_quantities(
    quantities=columns, 
    filters=galaxy_cuts, 
    native_filters=tract_filter([3830])#, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
)

extragalactic_cut = extragalactic_cat.get_quantities(
    quantities=columns_truth,
    filters=cosmoDC2_cuts,
    ) 

#Data_star_cut = Truth_match_cat.Get_quantities(
#    Quantities=Columns, 
#    Filters=Star_cuts, 
#    Native_filters=Tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#Data_extended_cut = Truth_match_cat.Get_quantities(
#    Quantities=Columns, 
#    Filters=Is_extended_cuts, 
#    Native_filters=Tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
#)

#sys.exit()

###Photoz plot
#fig = plt.figure(figsize=(12,8))
#plt.hist(x, range = (0, 3), bins = 30)
plt.hist(data_galaxy_cut['photoz_mode'], range = (0, 3), bins = 30, label="photo-z mode", histtype='step', color = 'blue');
#plt.hist(extragalactic_cut['redshift'], range = (0, 3), bins = 30, label="photo-z true", histtype='step', color = 'red');
#plt.plot(cat.photoz_pdf_bin_centers, sumpdf*3.,label="summed $p(z)$",lw=2,c='r');
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
plt.xlabel("redshift");
#plt.legend([g1,g2], ["dr3", "cosmoDC2"], title = 'Photo-z', loc='upper right') 
plt.legend(title = 'Photo-z', loc='upper right') 
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/photoz/"
#sys.mkdir(outpath);
plt.savefig(outpath+"photoz_mode.png", bbox_inches='tight')
plt.close()
sys.exit()











###Compute sky area
d_ra = data["ra"].max() - data["ra"].min()
d_dec = data["dec"].max() - data["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(data["dec"])))
sky_area_sq_deg = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg)
print('********Plot********')
###Data plot
mag_bins = np.linspace(14, 24, 200)
cdf = np.searchsorted(data["mag_i_cModel"], mag_bins, sorter=data["mag_i_cModel"].argsort())
cdf_galaxy_cut = np.searchsorted(data_galaxy_cut["mag_i_cModel"], mag_bins, sorter=data_galaxy_cut["mag_i_cModel"].argsort())
cdf_star_cut = np.searchsorted(data_star_cut["mag_i_cModel"], mag_bins, sorter=data_star_cut["mag_i_cModel"].argsort())
cdf_extended_cut = np.searchsorted(data_extended_cut["mag_i_cModel"], mag_bins, sorter=data_extended_cut["mag_i_cModel"].argsort())
g1, = plt.semilogy(mag_bins, cdf / sky_area_sq_deg, color = 'black')
g2, = plt.semilogy(mag_bins, cdf_galaxy_cut / sky_area_sq_deg, color = 'red')
#g3, = plt.semilogy(mag_bins, cdf_star_cut / sky_area_sq_deg, color = 'blue')
#g4, = plt.semilogy(mag_bins, cdf_extended_cut / sky_area_sq_deg, color = 'orange')
###cosmoDC2
Cdf = np.searchsorted(truth["mag_i"], mag_bins, sorter=truth["mag_i"].argsort())
###Compute sky area cosmoDC2
d_ra = truth["ra"].max() - truth["ra"].min()
d_dec = truth["dec"].max() - truth["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(truth["dec"])))
sky_area_sq_deg_cosmoDC2 = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg_cosmoDC2)
g5, = plt.semilogy(mag_bins, Cdf / sky_area_sq_deg_cosmoDC2, color = 'blue')
plt.xlabel("$i$-band mag");
plt.ylabel("Galaxies /deg2 /0.05mag");
#plt.legend([g1,g2,g3,g4],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type=2)", "dr6 (ext=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
#plt.legend([g1,g2,g3],["dr6", "dr6 (truth_type=1)", "dr6 (truth_type!=1)"], title = 'DC2 (9 tracts) / is_good_match', loc='upper left')
plt.legend([g1,g2,g5],["dr6", "dr6 (ext=1)", "cosmoDC2"], title = 'DC2 (9 tracts)', loc='upper left')
#plt.grid();
#plt.grid(True, which="both", ls="-", )
plt.grid(which='major', axis='both', linestyle='-', linewidth='0.5', color='grey')
# Customize the minor grid
plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.5', color='grey')
###Define the output path for figures
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/matching/"
#sys.mkdir(outpath);
plt.savefig(outpath+"dr6.png", bbox_inches='tight') 
plt.close()
print('********Plot saved********')

###CosmoDC2 plot
# Cdf = Np.Searchsorted(Truth["Mag_i_cModel"], Mag_bins, Sorter=Truth["Mag_i_cModel"].Argsort())
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
