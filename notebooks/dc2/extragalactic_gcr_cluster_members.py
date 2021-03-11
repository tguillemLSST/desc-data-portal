#!/usr/bin/env python
# coding: utf-8

# # Plotting Positions of Galaxy Cluster Members in Extragalactic Catalogs
# 
# In this example script we show how to select for and look at members of individual clusters.
# 
# Owners: **Dan Korytov [@dkorytov](https://github.com/LSSTDESC/DC2-analysis/issues/new?body=@dkorytov)**, **Patricia Larsen**
# 
# Last verified run: **Nov 30, 2018** by @yymao
# 
# This notebook demonstrates how to access the extra galactic catalog through the Generic Catalog Reader (GCR, https://github.com/yymao/generic-catalog-reader) as well as how filter on galaxy features and cluster membership.
# 
# __Objectives__:
# 
# After working through and studying this Notebook you should be able to
# 
# 1. Access extragalactic catalogs (protoDC2, cosmoDC2) through the GCR.
# 2. Select galaxy cluster centrals as a proxy for clusters.
# 3. Select galaxies in individual clusters by using the host_id quantity.
# 4. Plotting galaxy clustermembers positions on the sky as well as their comoving position in space.
# 
# 
# __Logistics__: This notebook is intended to be run through the JupyterHub NERSC interface available here: https://jupyter-dev.nersc.gov. To setup your NERSC environment, please follow the instructions available here: https://confluence.slac.stanford.edu/display/LSSTDESC/Using+Jupyter-dev+at+NERSC

# In[ ]:


import GCRCatalogs
import numpy as np
from astropy.table import Table
from GCR import GCRQuery
import matplotlib.pyplot as plt
import sys

# ### Reading catalog
# We load in the catalog with the "load_catalog" command, and then the values with the "get_quantities" command using filters to select sub-samples of the catalog.  
# 
# ### Help for error messages:
# If this fails to find the appropriate quantities, check that the desc-python kernel is being used and if this is not available source the kernels by running the following command on a terminal at nersc: "source /global/common/software/lsst/common/miniconda/setup_current_python.sh"
# 
# We are loading in a smaller version of the full cosmoDC2 catalog - this contains the same information as the full catalog but with a smaller sky area.

# In[ ]:


gc = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_small')
print('extragalactic_cat quantities')
print(sorted(gc.list_all_quantities()))
#for qty in ['cosmodc2_hp', 'cosmodc2_id', 'dec', 'flux_g', 'flux_g_noMW', 'flux_i', 'flux_i_noMW', 'flux_r', 'flux_r_noMW', 'flux_u', 'flux_u_noMW', 'flux_y', 'flux_y_noMW', 'flux_z', 'flux_z_noMW', 'host_galaxy', 'id', 'is_pointsource', 'is_variable', 'mag_g', 'mag_g_noMW', 'mag_i', 'mag_i_noMW', 'mag_r', 'mag_r_noMW', 'mag_u', 'mag_u_noMW', 'mag_y', 'mag_y_noMW', 'mag_z', 'mag_z_noMW', 'patch', 'ra', 'redshift', 'tract', 'truth_type', 'redshift', 'redshift_true']:
for qty in ['redshift', 'redshift_true']:
    info_dict = gc.get_quantity_info(qty)
    print(qty,info_dict)

###Galaxy and cluster selection
galaxy_data = gc.get_quantities(['ra', 'dec', 'mag_g', 'mag_r', 'mag_i', 'halo_id', 'redshift'], filters=['mag_i < 24'])
cluster_data = gc.get_quantities(['ra','dec', 'halo_mass', 'halo_id', 'redshift'], 
                                 filters=['is_central', 'halo_mass > 1e14', 'redshift < 0.3'])

###Define the output path for figures
outpath = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/cosmoDC2/"

#2D plot mass VS redshift for clusters
plt.figure()
plt.scatter(
    cluster_data['redshift'], 
    cluster_data['halo_mass'],
    marker='.',
    c='black'#,
    #s='5'
    )
plt.semilogy()
plt.legend(loc='best', framealpha=0.3)
plt.xlabel('z')
plt.ylabel('M ')
#plt.title('Halo ID:  {}\nHalo Mass:  {:.2e} h^-1 Msun'.format(cluster['halo_id'], cluster['halo_mass']))
plt.savefig(outpath+"mass_redshift.png", bbox_inches='tight')
plt.close()
print('********Plot saved********')
  
#number of clusters to debug
nmax = 3

cluster_data = Table(cluster_data)
for i, cluster in enumerate(cluster_data):
    if (i >= nmax):
        break # plot only the first 3
    members = GCRQuery('halo_id == {}'.format(cluster['halo_id'])).filter(galaxy_data)
    plt.figure()
    plt.scatter(
        members['ra'], 
        members['dec'], 
        s=(19-members['mag_r'])*8, 
        label='Galaxy Members [{}]'.format(len(members['ra']))
    )
    plt.plot(cluster['ra'], cluster['dec'], 'xr', label='Cluster Center')
    plt.legend(loc='best', framealpha=0.3)
    plt.xlabel(r'ra [deg]')
    plt.ylabel(r'dec [deg]')
    plt.title('Halo ID:  {}\nHalo Mass:  {:.2e} h^-1 Msun'.format(cluster['halo_id'], cluster['halo_mass']))
    plt.savefig(outpath+format(cluster['halo_id'])+".png", bbox_inches='tight')
    plt.close()
    print('********Plot saved********')
        
#    plt.show()

###Color diagrams
outpath_colors = "/sps/lsst/users/tguillem/DESC/config/tests_DC2/desc-data-portal/notebooks/dc2/plots/cosmoDC2/colors/"
for i, cluster in enumerate(cluster_data):
    if (i >= nmax):
        break # plot only the first 3
    members = GCRQuery('halo_id == {}'.format(cluster['halo_id'])).filter(galaxy_data)
    #plot g-r vs r-i
    plt.figure()
    plt.scatter(members['mag_g'] - members['mag_r'],
                members['mag_r'] - members['mag_i'],
                marker='.',
                label='galaxies')
    plt.xlabel('g-r')
    plt.ylabel('r-i')
    #plt.xlim(-1, +2)
    #plt.ylim(25, 15)
    plt.title('')
    plt.savefig(outpath_colors+format(cluster['halo_id'])+".png", bbox_inches='tight')
    plt.close()

    #
    #plot r-i vs i
    plt.figure()
    plt.scatter(members['mag_i'],
                members['mag_r']-members['mag_i'],
                marker='.',
                label='galaxies')
    plt.xlabel('i')
    plt.ylabel('r-i')
    #plt.xlim(-1, +2)
    #plt.ylim(25, 15)
    plt.title('')
    plt.savefig(outpath_colors+format(cluster['halo_id'])+"_rel.png", bbox_inches='tight')
    plt.close()
    print('********Plot saved********')

sys.exit()
#########################################################################################
# ### Extensions:
# We can load further information on the cluster members. For instance the second cluster looks a little odd in projected space, so we re-make these plots in comoving cartesian coordinates x and y in the example below. We also map the colours to the x-direction velocities. 
# 
# To do this you need to load the required quantities from the catalog before using them. A simple way to double check the quantity names is the command "gc.list_all_quantities()".

# In[ ]:


galaxy_data = gc.get_quantities(['ra', 'dec', 'mag_r', 'halo_id', 'position_x', 'position_y', 'velocity_x', 'velocity_y'], filters=['mag_r < 19'])
cluster_data = gc.get_quantities(['ra','dec', 'halo_mass', 'halo_id', 'position_x', 'position_y'], 
                                 filters=['is_central', 'halo_mass > 1e14', 'redshift < 0.2'])


# In[ ]:


cluster_data = Table(cluster_data)
for i, cluster in enumerate(cluster_data):
    if (i >= 3):
        break # plot only the first 3
    members = GCRQuery('halo_id == {}'.format(cluster['halo_id'])).filter(galaxy_data)
    plt.figure()
    plt.scatter(
        members['position_x'],
        members['position_y'],
        s=(19-members['mag_r'])*8, 
        label='Galaxy Members [{}]'.format(len(members['ra'])),
        c=members['velocity_x'],
        cmap='viridis')
    plt.plot(cluster['position_x'], cluster['position_y'],'xr',label='Cluster Center', alpha=0.8)
    plt.legend(loc='best',framealpha=0.3)
    plt.xlabel(r'x [Mpc/h]')
    plt.ylabel(r'y [Mpc/h]')
    plt.title('Halo ID: {}\nHalo Mass: {:.2e}'.format(cluster['halo_id'], cluster['halo_mass']))
    plt.colorbar(label='km/s')
plt.show()

