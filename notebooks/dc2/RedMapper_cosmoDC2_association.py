#!/usr/bin/env python
# coding: utf-8

### Adaptation of https://github.com/LSSTDESC/CLUSTER_VALIDATION [M. Ricci]

###import
import GCRCatalogs
from GCRCatalogs.helpers.tract_catalogs import tract_filter, sample_filter

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import sys
import os

from cluster_validation.opening_catalogs_functions import *
from cluster_validation.association_methods import *
from cluster_validation.plotting_functions import *
from cluster_validation.association_statistics import *

plt.rcParams['figure.figsize'] = [9.5, 6]
plt.rcParams.update({'font.size': 18})
#plt.rcParams['figure.figsize'] = [10, 8] for big figures

outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/redmapper/"
if not os.path.exists(outpath):
    os.mkdir(outpath)

####function to open truth and detection catalogs
RM_cat_name = 'cosmoDC2_v1.1.4_redmapper_v0.5.7'
DC2_cat_name = 'cosmoDC2_v1.1.4'
#DC2_cat_name = 'cosmoDC2_v1.1.4_small'       

min_richness = 20
min_halo_mass = 3e14 #Msun

cluster_data, member_data, truth_data, gc, gc_truth = RM_DC2_cat_open(RM_cat_name,DC2_cat_name,min_richness, min_halo_mass, cluster_only=False)

#take only the halo
halo_data = truth_data[truth_data['is_central']==True]

min_richness = 25
    # Get the redMaPPer catalog
gc = GCRCatalogs.load_catalog(RM_cat_name)
    # Select out the cluster and member quantities into different lists
quantities = gc.list_all_quantities()
cluster_quantities = [q for q in quantities if 'member' not in q]
member_quantities = [q for q in quantities if 'member' in q]
    
# Read in the cluster and member data
#reduce to _small
ra_min, ra_max = 60, 70
dec_min, dec_max = -46, -34
query = GCRCatalogs.GCRQuery('richness > ' + str(min_richness))
                             #'ra > {}'.format(ra_min),
                             #'ra < {}'.format(ra_max),
                             #'dec > {}'.format(dec_min),
                             #'dec < {}'.format(dec_max)
                             #)
cluster_data = Table(gc.get_quantities(cluster_quantities, [query]))
member_data = Table(gc.get_quantities(member_quantities))

####General catalog properties
print("Number of elements in the truth catalog = ", len(truth_data))
print("Number of halos in the truth catalog = ", len(halo_data))
print("Number of clusters in the detection catalog = ", len(cluster_data))
print("Cluster catalog sky area = ", gc.sky_area, "deg2")
print("Truth catalog sky area = ", gc_truth.sky_area, "deg2")

###Define same cosmological parameters as in the truth catalog (cosmoDC2)
cosmo = gc_truth.cosmology
#check the cosmological parameters in the two catalogs
#print('Cosmo in truth catalog:', gc_truth.cosmology)
#print('Cosmo in detection catalog:', gc_truth.cosmology)

###Basic visualization
plt.figure()
plt.plot(cluster_data['ra_cen_0'],cluster_data['dec_cen_0'],'b.')
plt.plot(truth_data['ra'],truth_data['dec'],'rx',alpha=0.1)
plt.xlabel("ra")
plt.ylabel("dec");
plt.savefig(outpath+"clusters.png", bbox_inches='tight')
plt.close()
print('********Plot saved********')
#sys.exit()

###Associate Redmapper detections to true DC2 halos

# **Association method scheme :**
# - search for RM associations in a given cylinder around each DC2 object (radius $\theta_{max}$ and width $\Delta_z\times(1+z_{object})$)  
# ($\theta_{max}$  can be **fixed** - in Mpc or arcmin - **or scale** with halo mass and richness and $\Delta_z$ can be infinite)
# - select the **nearest** match in projected distance **or** the one with the **more galaxies in common**
# - repeat in the other direction (RM>DC2)
# - return the associations which are bijective
# 
# 

### - 1rst association criteria : nearest within a cylinder with $\Delta_z=\inf, \theta_{max} = 1$ arcmin
#criteria
delta_zmax = np.inf
theta_max = 1. #arcmin
theta_max_type = "fixed_angle"
method = "nearest"

match_num_1w, match_num_2w, ind_bij = volume_match(halo_data, cluster_data, delta_zmax, theta_max, theta_max_type, method, cosmo, truth_data, member_data)
#truth_to_det_match_numbers, det_to_truth_match_number, bijective_match_indices

#statistics
print ("number of bijective associations", number_of_associations(ind_bij))
print ("number and fraction of fragmentation", fragmentation(match_num_1w, ind_bij, method="bij"))
print ("number and fraction of overmerging", overmerging(match_num_2w, ind_bij, method="bij"))
print ("completeness", completeness(halo_data, ind_bij, gc, gc_truth))
print ("purity", purity(cluster_data, ind_bij, gc, gc_truth))
#sys.exit()

#check plot
f,a,b = plot_cluster_and_halo_position(halo_data, cluster_data, match_num_1w, match_num_2w, ind_bij, outpath)


# In[ ]:


#plot to check redshift
plot_redshift_comparison(halo_data, cluster_data, ind_bij, outpath)
sys.exit()

# ### - 2nd association criteria : nearest within a cylinder with $\Delta_z=0.05, \theta_{max} = 1 Mpc$ 

# In[14]:


#criteria
delta_zmax = 0.05
theta_max = 1. #Mpc
theta_max_type = "fixed_dist"
method = "nearest"

match_num_1w, match_num_2w, ind_bij = volume_match(halo_data, cluster_data, delta_zmax, theta_max, theta_max_type, method, cosmo, truth_data, member_data)
#truth_to_det_match_numbers, det_to_truth_match_number, bijective_match_indices


# In[15]:


#statistics
print ("number of bijective associations", number_of_associations(ind_bij))
print ("number and fraction of fragmentation", fragmentation(match_num_1w, ind_bij, method="bij"))
print ("number and fraction of overmerging", overmerging(match_num_2w, ind_bij, method="bij"))
print ("completeness", completeness(halo_data, ind_bij, gc, gc_truth))
print ("purity", purity(cluster_data, ind_bij, gc, gc_truth))


# In[ ]:


#check plot
f,a,b = plot_cluster_and_halo_position(halo_data, cluster_data, match_num_1w, match_num_2w, ind_bij)


# ### - 3rd association criteria : nearest within a cylinder with $\Delta_z=0.05, \theta_{max} = R(M_{halo}, \lambda_{cluster})$ 

# In[16]:


#criteria
delta_zmax = 0.05
theta_max = 1. #Mpc
theta_max_type = "scaled"
method = "nearest"

match_num_1w, match_num_2w, ind_bij = volume_match(halo_data, cluster_data, delta_zmax, theta_max, theta_max_type, method, cosmo, truth_data, member_data)
#truth_to_det_match_numbers, det_to_truth_match_number, bijective_match_indices


# In[17]:


#statistics
print ("number of bijective associations", number_of_associations(ind_bij))
print ("number and fraction of fragmentation", fragmentation(match_num_1w, ind_bij, method="bij"))
print ("number and fraction of overmerging", overmerging(match_num_2w, ind_bij, method="bij"))
print ("completeness", completeness(halo_data, ind_bij, gc, gc_truth))
print ("purity", purity(cluster_data, ind_bij, gc, gc_truth))


# In[18]:


#check plot
f,a,b = plot_cluster_and_halo_position(halo_data, cluster_data, match_num_1w, match_num_2w, ind_bij)


# ### - 4rth association criteria : highest common membership within a cylinder with $\Delta_z=0.05, \theta_{max} = 1 Mpc$ 

# In[19]:


#criteria
delta_zmax = 0.05
theta_max = 1. #Mpc
theta_max_type = "fixed_dist"
method = "membership"

match_num_1w, match_num_2w, ind_bij = volume_match(halo_data, cluster_data, delta_zmax, theta_max, theta_max_type, method, cosmo, truth_data, member_data)
#truth_to_det_match_numbers, det_to_truth_match_number, bijective_match_indices


# In[20]:


#statistics
print ("number of bijective associations", number_of_associations(ind_bij))
print ("number and fraction of fragmentation", fragmentation(match_num_1w, ind_bij, method="bij"))
print ("number and fraction of overmerging", overmerging(match_num_2w, ind_bij, method="bij"))
print ("completeness", completeness(halo_data, ind_bij, gc, gc_truth))
print ("purity", purity(cluster_data, ind_bij, gc, gc_truth))


# In[21]:


#check plot
f,a,b = plot_cluster_and_halo_position(halo_data, cluster_data, match_num_1w, match_num_2w, ind_bij)


# In[ ]:





# In[ ]:




