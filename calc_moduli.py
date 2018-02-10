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
import scipy.io

k_b = 1.38064852e-23 # Boltzmann constant	
T=313
sumcoll =0;
gridcoll =0;
test_frame = frame(2*4096)
number_of_frames = 1500
number_of_blocks = 1
frames_per_block = int(number_of_frames/number_of_blocks)
M=62
dir_all = np.zeros((M,M,3,number_of_blocks))
n_long=[]
n_trans=[]
s_wavenum = []
for block in range(number_of_blocks):    
    sumcoll =0;
    gridcoll =0;
    for i in range(frames_per_block):
        test_frame.load_next_frame()
       
        lipid_grid = gf.create_lipid_grid(test_frame.bilayer.upper,M, 
                                      test_frame.box_size)
        director_grid = gf.create_director_grid(lipid_grid,M, 
                                                test_frame.box_size)
        lipid_grid2 = gf.create_lipid_grid(test_frame.bilayer.lower,M, 
                                      test_frame.box_size)
        director_grid2 = gf.create_director_grid(lipid_grid2,M, 
                                                 test_frame.box_size)
        director_grid3 = 0.5*(director_grid-1*director_grid2)
        gf.normalize_grid(director_grid3,M)
        q_grid = gf.create_grid_qvalues(M,test_frame.box_size)
        n_q=gf.fourier_trans_grid(director_grid3,M,test_frame.box_size[0] )
        n_pow_sum, w_grid_sum = gf.collect(n_q,q_grid,M,test_frame.box_size[0])
        sumcoll += n_pow_sum/frames_per_block
        gridcoll += w_grid_sum/frames_per_block
        sorted_wavenum_u, n_long_u, n_trans_u = (
                gf.get_moduli(gridcoll,sumcoll[:,:,0],sumcoll[:,:,1]))
        print ("\r" + "{0:.5f}".format(100*(i/frames_per_block)) + "% done \r",end="")
    dir_all[:,:,:,block] = director_grid3
    n_long.append(n_long_u)
    n_trans.append(n_trans_u)
    s_wavenum.append(sorted_wavenum_u)
    
        
        
    print ("\r" + "{0:.5f}".format(100*(block/number_of_blocks)) + "% done \r",end="")



#scipy.io.savemat('dir_all_filt_Long.mat', mdict={'dir_all': dir_all})
K_C_all = np.zeros(number_of_blocks)
for i in range(number_of_blocks):
    K_C_tmp = np.square(s_wavenum[i])*n_trans[i]
    K_C_all[i] = k_b*T/np.mean(K_C_tmp[0:2])

sorted_wavenum_u, n_long_u, n_trans_u = gf.get_moduli(gridcoll,sumcoll[:,:,0],sumcoll[:,:,1])
plt.plot(sorted_wavenum_u,n_long_u)
K_C = np.square(sorted_wavenum_u)*n_trans_u
K_C = k_b*T/np.mean(K_C[0:2])
print('Kc = ',K_C)
plt.show()