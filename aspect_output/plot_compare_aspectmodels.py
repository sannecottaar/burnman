#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    # Reads aspect output and/or computed seismic velocities for up to four models and plots these
    # Aspect output file and h5 file with computed seismic velocities (computed using compute_velocities_aspectoutput.py) are expected in output/#model_name/solution/
    # Figures are saved as 'all_models_'+value_to_plot+'.pdf'
"""




import sys, os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import average_solid_melt

# hack to allow scripts to be placed in subdirectories next to burnman:
sys.path.insert(1, os.path.abspath('../'))
sys.path.insert(1, os.path.abspath('../boukare/'))

# Here we import the relevant modules from BurnMan.
import burnman
from burnman import minerals
from burnman import averaging_schemes
import ulvz_melt_model
import glob

'''
input
'''
# Models to plot
model_names = ['reference','no_motion_of_melt', 'melt_is_less_dense','with_heterogeneities']

# Plot value from aspect output across model
plot_aspect_output = False
value_to_plot = 'temperature' # implemented are 'temperature', 'melt', 'solid_fe'

# Plot computed Vs across models
plot_seismic_velocity = True


'''
'''

plt.figure(figsize=(10,6))


if plot_aspect_output:
    for m,model_name in enumerate(model_names):
        # Define mesh and solution file
        mesh_file = glob.glob('output/'+model_name+'/solution/mesh-*.h5')[0]
        solution_file = glob.glob('output/'+model_name+'/solution/solution-*.h5')[0]

        ## Load the mesh
        mesh = h5.File(mesh_file,'r')
        nodes = np.array(mesh['nodes'])
        
        ## Load aspct output
        solution = h5.File(solution_file,'r')
        print('keys in solution', list(solution.keys()))


        ## identify duplicate nodes - note that this is a bit sloppy because it presumes that the nodal values are also duplicated, which may not be true for DG fields
        unique_nodes, unique_indices, unique_inverse, unique_counts = np.unique(nodes,axis=0,return_index=True,return_inverse=True,return_counts=True)

        x = unique_nodes[:,0]
        y = unique_nodes[:,1]

        # load T and P
        temperatures = np.array(solution['T'])[unique_indices][:,0]
        pressures = np.array(solution['p_f'])[unique_indices][:,0] # 'p' goes negative, 'p_f' does not

        #viscosity and density
        viscosity = np.array(solution['viscosity'])[unique_indices][:,0]
        density = np.array(solution['density'])[unique_indices][:,0]

        # load melt fraction
        porosity  = np.array(solution['porosity'])[unique_indices][:,0] # porosity can go slightly negative
        melt_frac = np.array(solution['melt_fraction'])[unique_indices][:,0] # runs from 0-1

        # load composition
        bulk_composition = np.array(solution['bulk_composition'])[unique_indices][:,0]
        melt_fe = np.array(solution['molar_Fe_in_melt'])[unique_indices][:,0]
        solid_fe = np.array(solution['molar_Fe_in_solid'])[unique_indices][:,0]

        if value_to_plot == 'temperature':
            # plot temperature
            plt.subplot(2,2,m+1)
            plt.tricontourf(x/1.e3,y/1.e3,temperatures, 100, cmap='hot', extend = 'both')
            plt.colorbar()
            #plt.tricontour(x/1.e3,y/1.e3,temperatures,levels=[ 0.92*3900,3900], colors='k')
            plt.ylabel('height above CMB (km)')
            plt.xlabel('(km)')
        
        if value_to_plot == 'melt':
            # plot melt
            plt.subplot(2,2,m+1)
            plot_val = np.linspace(0,.1, 21, endpoint=True)
            pl = plt.tricontourf(x/1.e3,y/1.e3, melt_frac , plot_val,  cmap='copper', vmin =0, vmax=0.1, extend = 'both')
            plt.colorbar()
            #plt.tricontour(x/1.e3,y/1.e3,temperatures,levels=[ 0.92*3900,3900], colors='k')
            plt.ylabel('height above CMB (km)')
            plt.xlabel('(km)')
        
        if value_to_plot == 'solid_fe':
            # plot iron content in nsolid
            plt.subplot(2,2,m+1)
            plot_val = np.linspace(0,.2, 21, endpoint=True)
            pl = plt.tricontourf(x/1.e3,y/1.e3, solid_fe , plot_val,  cmap='copper_r', vmin =0, vmax=0.2, extend = 'both')
            plt.colorbar()
            #plt.tricontour(x/1.e3,y/1.e3,temperatures,levels=[ 0.92*3900,3900], colors='k')
            plt.ylabel('height above CMB (km)')
            plt.xlabel('(km)')
        
        
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(model_name)
    
    plt.savefig('all_models_'+value_to_plot+'.pdf')
    plt.show()


if plot_seismic_velocity:
    ## Load and plot seismic velocities
    plt.figure(figsize=(10,6))


        vel_file = glob.glob('output/'+model_name+'/solution/seismic_velocities.h5')[0]
        vels = h5.File(vel_file,'r')
        vs = np.array(vels['vs'])
        x= np.array(vels['x'])
        y = np.array(vels['y'])
        dvs = (np.array(vs)/np.max(vs)-1.)*100.

        plot_val = np.linspace(-20,0, 21, endpoint=True)

        # Plot velocity
        plt.subplot(2,2,m+1)
        pl= plt.tricontourf(x/1.e3,y/1.e3,dvs, plot_val, cmap='OrRd_r', vmin = -20., vmax = 0, extend = 'both')
        plt.ylabel('height above CMB (km)')
        plt.xlabel('(km)')
        bar= plt.colorbar()
        pl.set_clim(-20., 0)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(model_name)

#


    plt.savefig("all_models_vs.pdf")
    plt.show()


####
