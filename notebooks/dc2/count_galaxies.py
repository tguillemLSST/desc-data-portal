#!/usr/bin/env python
# coding: utf-8

# # Rubin LSST DESC DC2: Accessing Object Table with GCRCatalogs
# 
# **Authors**: Yao-Yuan Mao (@yymao), Francois Lanusse (@EiffL), Javier Sanchez (@fjaviersanchez), Michael Wood-Vasey (@wmwv), Rachel Mandelbaum (@rmandelb)
# 
# This notebook will illustrate the basics of accessing the Object Table, which contains the detected objects at the coadd level using GCRCatalogs.
# 
# **Learning objectives**: After going through this notebook, you should be able to:
#   1. Load and efficiently access a DC2 object table with the GCRCatalogs
#   2. Understand and have references for the object table schema
#   3. Apply cuts to the catalog using GCRQuery
#   4. Have an example of quality cuts and simple star/galaxy separation cut

# ## Before you start
# 
# Make sure you have followed the instructions on the [DESC Data Portal](https://lsstdesc-portal.nersc.gov/) to 
# download the data files, install `GCRCatalogs`, and set up `root_dir` for `GCRCatalogs`.
# 
# In this example notebook we will use up to 4 tracts: 3830, 3831, 4028, 4029. 
# The example set from the Data Portal includes these four tracts, so If you have downloaded it you are all set. 
# If you have downloaded different tracts, remember to change the tract number(s) used in `tract_filter` when you run this notebook.

# ## Import necessary packages

# In[1]:


#empty


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
plt.ioff() 
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter
from GCRCatalogs import GCRQuery


# ## Access object table with GCRCatalogs
# 
# The [GCRCatalogs](https://github.com/LSSTDESC/gcr-catalogs) package is a DESC project which aims at gathering in one convenient location various simulation/data catalogs made available to the collaboration. In this section, we illustrate how to use this tool to access the object catalogs from DC2 Run2.2i.

# In[4]:


GCRCatalogs.get_public_catalog_names()


# In[5]:


obj_cat = GCRCatalogs.load_catalog("dc2_object_run2.2i_dr6_v2")


# ### Object Table schema
# 
# To see the "columns" (sometimes called "quantities") available in the object table, you can call `list_all_quantities()`.
# The returned list of `list_all_quantities` is not sorted. Sorting it would make it easier to read.

# In[6]:


sorted(obj_cat.list_all_quantities())


# The definitions, units, and types of these fields are documented in the data release note (Table 1).
# As explained in the note, the values exposed here are not the native columns produced by the LSST Science Pipelines.
# Instead, this schema strives to mimic the schema of the LSST Data Products Definition Document [DPDD](http://ls.st/dpdd), 
# which is an effort made by the Rubin Observatory to standardize the format of the official Data Release Products (DRP).
# 
# If the data release note is many clicks away, you can also check the definitions, units, and types of these columns
# using `get_quantity_info()`. Here's an example:

# In[7]:


#empty
extragalactic_cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small')
#extragalactic_cat = GCRCatalogs.load_catalog('skysim5000_v1.1.1_small')
sorted(obj_cat.list_all_quantities())

# In[10]:


#empty
#is_galaxy = GCRQuery('truth_type == 1')  # Extended objects (primarily galaxies)

#clean = GCRQuery('clean')  # The source has no flagged pixels (interpolated, saturated, edge, clipped...) 
#                           # and was not skipped by the deblender
#bright =  GCRQuery('mag_i>23')

cosmoDC2_cuts = [
    #is_galaxy,
    #bright,
    GCRQuery((np.isfinite, 'mag_i')),  # Select objects that have i-band cmodel magnitudes
]


# In[ ]:


###truth selection
columns = ["ra", "dec", "mag_i"]
truth = extragalactic_cat.get_quantities(
    quantities=columns,
    filters=cosmoDC2_cuts
)

###data selection
is_extended = GCRQuery('extendedness == 1')  # Extended objects (primarily galaxies)
clean = GCRQuery('clean')  # The source has no flagged pixels (interpolated, saturated, edge, clipped...) 
                           # and was not skipped by the deblender

galaxy_cuts = [
    clean,
    is_extended, 
    GCRQuery((np.isfinite, 'mag_i')),  # Select objects that have i-band cmodel magnitudes
]

# for purposes other than the cuts!
columns = ["ra", "dec", "mag_i"]

data = obj_cat.get_quantities(
    quantities=columns, 
    filters=galaxy_cuts, 
    native_filters=tract_filter([3830, 3831, 3832, 4028, 4029, 4030, 3636, 3637, 3638])
)

#compute sky are
d_ra = data["ra"].max() - data["ra"].min()
d_dec = data["dec"].max() - data["dec"].min()
cos_dec = np.cos(np.deg2rad(np.median(data["dec"])))
sky_area_sq_deg = (d_ra * cos_dec) * d_dec
print(sky_area_sq_deg)
print('********Plot********')
#data plot
mag_bins = np.linspace(14, 24, 200)
cdf = np.searchsorted(data["mag_i"], mag_bins, sorter=data["mag_i"].argsort())
#plt.hist(values, num_bins)
plt.semilogy(mag_bins, cdf / sky_area_sq_deg)
plt.xlabel("$i$-band mag");
plt.ylabel("Galaxies /deg2 /0.05mag");

plt.grid();  # add grid to guide our eyes

#Define the output path for figures
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/"
plt.savefig(outpath+"data.png", bbox_inches='tight') 

#cosmoDC2 plot
#cdf = np.searchsorted(truth["mag_i"], mag_bins, sorter=truth["mag_i"].argsort())
#sky_area_sq_deg = 50;
#plt.semilogy(mag_bins, cdf / sky_area_sq_deg)
#plt.xlabel("$i$-band mag");
#plt.ylabel("");

#plt.grid();  # add grid to guide our eyes

#Define the output path for figures
#outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/"
#plt.savefig(outpath+"truth", bbox_inches='tight') 
