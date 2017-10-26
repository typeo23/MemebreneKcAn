#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 20:44:21 2017

@author: itay
"""

import matplotlib.pyplot as plt
from frame import Frame as frame
import grid_functions as gf
import numpy as np

k_b = 1.38064852e-23 # Boltzmann constant	
T=313
sumcoll =0;
gridcoll =0;
test_frame = frame(1024)
nimber_of_frmes = 1
M=36
for i in range(nimber_of_frmes):
    test_frame.load_next_frame()
   
    lipid_grid = gf.create_lipid_grid_closest2(test_frame.bilayer.upper,M, 
                                  test_frame.box_size)
    director_grid = gf.create_director_grid(lipid_grid,M, 
                                            test_frame.box_size)
    lipid_grid2 = gf.create_lipid_grid_closest2(test_frame.bilayer.lower,M, 
                                  test_frame.box_size)
    director_grid2 = gf.create_director_grid(lipid_grid2,M, 
                                             test_frame.box_size)
    director_grid3 = 0.5*(director_grid-director_grid2)
    gf.normalize_grid(director_grid3,M)
    q_grid = gf.create_grid_qvalues(M,test_frame.box_size)
    n_q=gf.fourier_trans_grid(director_grid3,M,test_frame.box_size[0] )
    n_pow_sum, w_grid_sum = gf.collect(n_q,q_grid,M,test_frame.box_size[0])
    sumcoll += n_pow_sum/nimber_of_frmes
    gridcoll += w_grid_sum/nimber_of_frmes


sorted_wavenum_u, n_long_u, n_trans_u = gf.get_moduli(gridcoll,sumcoll[:,:,0],sumcoll[:,:,1])
plt.plot(sorted_wavenum_u,sorted_wavenum_u*sorted_wavenum_u*n_trans_u)
K_C = np.square(sorted_wavenum_u)*n_trans_u
K_C = k_b*T/np.mean(K_C[0:2])
print('Kc = ',K_C)
plt.show()