# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 20:36:33 2021

@author: 18313
"""
from aurora.solps import solps_case
from SOLPS_Scripts.TOOLS import SET_WDIR
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

BASE,TOP = SET_WDIR('','')

fig10,ax10 = plt.subplots(nrows=1, ncols=1)

ax10.set_xlabel(r'$\rho_p$')
ax10.set_ylabel(r'$n_n$')

'''
ax1[0].set_xlabel(r'$\rho_p$')
ax1[0].set_ylabel(r'$n_n$')

ax1[1].set_xlabel(r'$\rho_p$')
ax1[1].set_ylabel(r'$n_n/n_e$')
'''

label_list = [
    'C-Mod 1.3MA',
    'C-Mod 1.0MA',
    'DIII-D 175060',
    'ITER']

path_list = [
    '{}cmod/012home/Attempt65N/Output'.format(BASE),
    '{}cmod/025home/Attempt48N/Output'.format(BASE),
    '{}d3d/Attempt86/Output'.format(BASE),
    '{}iter/iter_solps_jl'.format(BASE)]

run_list = ['','','','orig_D1.95e23_Ne2.00e20.done.ids']

gfilepath_list = [
    '{}gfileProcessing/cmod_files/g1160718012.01209_974'.format(TOP),
    '{}gfileProcessing/cmod_files/g1160718025.01209_974'.format(TOP),
    '{}gfileProcessing/d3d_files/g175060.02512'.format(TOP),
    '{}iter/iter_solps_jl/baserun/gfile_iter'.format(BASE),]

rhop_fsa_min=0.8
rhop_fsa_max=1.5

for i in range(len(label_list)):
    
    so = solps_case(path_list[i], gfilepath_list[i], solps_run=run_list[i], form='full')

    # load neutral density
    
    rhop_fsa, nn_fsa, rhop_LFS, nn_LFS, rhop_HFS, nn_HFS = so.get_radial_prof(so.quants['nn'], dz_mm=1, plot=False)
    
    rhop_fsa, ne_fsa, rhop_LFS, ne_LFS, rhop_HFS, ne_HFS = so.get_radial_prof(so.quants['ne'], dz_mm=1, plot=False)
    
    ax10.semilogy(rhop_fsa, nn_fsa, label=label_list[i])
    #ax1[1].semilogy(rhop_LFS, nn_LFS/ne_LFS, label=label_list[i])
    
    if rhop_fsa[0]>rhop_fsa_min:
        rhop_fsa_min=rhop_fsa[0]
        
    if rhop_fsa[-1]<rhop_fsa_max:
        rhop_fsa_max=rhop_fsa[-1]

ax10.set_xlim(xmin=rhop_fsa_min,xmax=rhop_fsa_max)
#ax1[1].set_xlim(xmin=rhop_LFS_min,xmax=rhop_LFS_max) 
    
ax10.legend()
#ax1[1].legend()