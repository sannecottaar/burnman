#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    # Reads in seismic velocity output for multiple aspect runs
    # I use this to explore how velocity changes as function of temperature, melt_frac, iron_fraction etc.
"""


import glob
import ulvz_melt_model
from burnman import averaging_schemes
from burnman import minerals
import burnman
import sys
import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import average_solid_melt

# hack to allow scripts to be placed in subdirectories next to burnman:
sys.path.insert(1, os.path.abspath('../'))
sys.path.insert(1, os.path.abspath('../boukare/'))

# Here we import the relevant modules from BurnMan.


# Models to plot
model_names = [
    'reference',
    'no_motion_of_melt',
    'melt_is_less_dense',
    'with_heterogeneities']
plt.figure(figsize=(10, 6))

resample = 10

for m, model_name in enumerate(model_names):

    # Define mesh and solution file
    mesh_file = glob.glob('output/' + model_name + '/solution/mesh-*.h5')[0]
    solution_file = glob.glob(
        'output/' +
        model_name +
        '/solution/solution-*.h5')[0]

    # Load the mesh
    mesh = h5.File(mesh_file, 'r')
    nodes = np.array(mesh['nodes'])

    # Load aspct output
    solution = h5.File(solution_file, 'r')
    print('keys in solution', list(solution.keys()))

    # identify duplicate nodes - note that this is a bit sloppy because it
    # presumes that the nodal values are also duplicated, which may not be
    # true for DG fields
    unique_nodes, unique_indices, unique_inverse, unique_counts = np.unique(
        nodes, axis=0, return_index=True, return_inverse=True, return_counts=True)

    x = unique_nodes[:, 0]
    y = unique_nodes[:, 1]

    # load T and P
    temperatures = np.array(solution['T'])[unique_indices][:, 0]
    # 'p' goes negative, 'p_f' does not
    pressures = np.array(solution['p_f'])[unique_indices][:, 0]

    #viscosity and density
    viscosity = np.array(solution['viscosity'])[unique_indices][:, 0]
    density = np.array(solution['density'])[unique_indices][:, 0]

    # load melt fraction
    # porosity can go slightly negative
    porosity = np.array(solution['porosity'])[unique_indices][:, 0]
    melt_frac = np.array(solution['melt_fraction'])[
        unique_indices][:, 0]  # runs from 0-1

    # load composition
    bulk_composition = np.array(solution['bulk_composition'])[
        unique_indices][:, 0]
    melt_fe = np.array(solution['molar_Fe_in_melt'])[unique_indices][:, 0]
    solid_fe = np.array(solution['molar_Fe_in_solid'])[unique_indices][:, 0]

    vel_file = glob.glob(
        'output/' +
        model_name +
        '/solution/seismic_velocities_solidonly_res10.h5')[0]
    vels = h5.File(vel_file, 'r')
    vs = np.array(vels['vs'])
    vp = np.array(vels['vp'])
    x = np.array(vels['x'])
    y = np.array(vels['y'])
    dvs = (np.array(vs) / np.max(vs) - 1.) * 100.
    dvp = (np.array(vp) / np.max(vp) - 1.) * 100.
    plt.xlim([0.05, 0.2])

    plt.subplot(2, 2, m + 1)
    plt.plot(temperatures[::resample], dvs / dvp, '.b',
             markersize=0.2, label='including_melt')

    vel_file = glob.glob(
        'output/' +
        model_name +
        '/solution/seismic_velocities_res10.h5')[0]
    vels = h5.File(vel_file, 'r')
    vs = np.array(vels['vs'])
    vp = np.array(vels['vp'])
    x = np.array(vels['x'])
    y = np.array(vels['y'])
    dvs = (np.array(vs) / np.max(vs) - 1.) * 100.
    dvp = (np.array(vp) / np.max(vp) - 1.) * 100.

    plt.subplot(2, 2, m + 1)
    plt.plot(temperatures[::resample], dvs / dvp, '.r',
             markersize=0.2, label='including_melt')
    plt.ylabel('Vs')
    plt.xlabel('melt_fraction')
    plt.xlim([0.0, 0.2])
    plt.ylim([1, 3])
    plt.title(model_names[m])

    if m == 2:
        plt.legend(loc=1)

plt.savefig("all_models_temperatures_dVsdVPratio.pdf")
plt.show()


####
