#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 01:27:00 2018

@author: itay

Run the main loop of the code and saves the output
"""

import matplotlib.pyplot as plt
from frame import Frame as frame
import grid_functions as gf
import numpy as np
import scipy.io
import pandas as pd

def clac_moduli(params):
    

    k_b = 1.38064852e-23 # Boltzmann constant	
    T=313 #temperatue
    
    #init the frmae with the approriate number of lipis
    test_frame = frame(params.lipids)
    number_of_frames = params.frames
    number_of_blocks = 1
    frames_per_block = int(number_of_frames/number_of_blocks)
    M=params.M
    
    # Accumulators
    dir_all = np.zeros((M,M,3,number_of_blocks)) #array conating all the direcotr grids
    n_long=[] #lonitual director componants spectra
    n_trans=[] #tansverse director component spectra
    n_long_com=[]
    n_trans_com =[]
    s_wavenum = [] #flattned wavenumbers
    
    # Begin the main loop on the number of blocks
    for block in range(number_of_blocks):    
        sumcoll =0;
        gridcoll =0;
        # Looping on all the  frames in a block
        for i in range(frames_per_block):
            #Load the next frame of lipids
            test_frame.load_next_frame()
            
            # save the size of the first box as q values will be caculated by this
            if (i==0):
                box_size = test_frame.box_size
           
            # producing wave number(q) values for the grid
            q_grid = gf.create_grid_qvalues(M,box_size)
            
            # Grdding the directors from the upper leadlet
            lipid_grid_upper = gf.create_lipid_grid(test_frame.bilayer.upper,M, 
                                          test_frame.box_size)
            # Director grid for the upper leaflet
            director_grid_upper = gf.create_director_grid(lipid_grid_upper,M, 
                                                    test_frame.box_size)
            #Gridding lipids for the lower leflet
            lipid_grid_lower = gf.create_lipid_grid(test_frame.bilayer.lower,M, 
                                          test_frame.box_size)
            #director grid for the lower leaflet
            director_grid_lower = gf.create_director_grid(lipid_grid_lower,M, 
                                                     test_frame.box_size)
            #Avragin the director fields 
            director_grid_av = 0.5*(director_grid_upper-director_grid_lower)
            
            #Nomalizing the directos field (to Z=1)
            gf.normalize_grid(director_grid_av,M)
            
            n_q=gf.fourier_trans_grid(director_grid_av,M,test_frame.box_size[0])
            #collecting the spectra
            n_pow_sum, w_grid_sum, n_comp = gf.collect(n_q,q_grid,M,box_size[0])
            sumcoll = n_pow_sum #/frames_per_block
            gridcoll = w_grid_sum #/frames_per_block
            # and a 1D array containing the q values and correspondant logitudenal/
            # Transverse componants of the director field
            sorted_wavenum_u, n_long_u, n_trans_u, n_long_c, n_trans_c = (
                    gf.get_moduli(gridcoll,sumcoll[:,:,0],sumcoll[:,:,1],
                                  n_comp[:,:,0], n_comp[:,:,1]))
            n_long.append(n_long_u)
            n_trans.append(n_trans_u)
            n_long_com.append(n_long_c)
            n_trans_com.append(n_trans_c)
            s_wavenum.append(sorted_wavenum_u)
            print ("\r" + "{0:.5f}".format(100*(i/frames_per_block)) 
                   + "% done \r",end="")
#        dir_all[:,:,:,block] = director_grid3
#        n_long.append(n_long_u)
#        n_trans.append(n_trans_u)
#        s_wavenum.append(sorted_wavenum_u)
        
            
            
        #print ("\r" + "{0:.5f}".format(100*(block/number_of_blocks)) + "% done \r",end="")
    
    
    
    scipy.io.savemat('dir_all_filt_Long.mat', mdict={'n_long': n_long,
                                                     'n_trans': n_trans,
                                                     's_wavenum': s_wavenum})
    #creating and saving pandas dataframe
    trans = np.array(n_trans)
    long = np.array(n_long)
    long_c = np.array(n_long_com)
    trans_c = np.array(n_trans_com)
    
    pd_trans = pd.DataFrame(trans, columns= s_wavenum[0])
    pd_long = pd.DataFrame(long, columns= s_wavenum[0])
    pd_trans.to_csv('trans.csv')
    pd_long.to_csv('long.csv')
    
    # Saving csv files
    np.savetxt('trans.csv',trans, delimiter = ',')
    np.savetxt('long.csv', long, delimiter=',')
    np.save('trans_c.csv',trans_c)
    np.save('long_c.csv', long_c)
    
    np.savetxt('wavenum.csv', s_wavenum[0], delimiter=',')
     
    
    K_C_all = np.zeros(number_of_blocks)
    for i in range(number_of_blocks):
        K_C_tmp = np.square(s_wavenum[i])*n_trans[i]
        K_C_all[i] = k_b*T/np.mean(K_C_tmp[0:2])
    
    sorted_wavenum_u, n_long_u, n_trans_u, nlc, ntc = gf.get_moduli(
                                            gridcoll,sumcoll[:,:,0],
                                            sumcoll[:,:,1],
                                            n_comp[:,:,0], n_comp[:,:,1])
    plt.plot(sorted_wavenum_u,n_long_u)
    K_C = np.square(sorted_wavenum_u)*n_trans_u
    K_C = k_b*T/np.mean(K_C[0:2])
    print('Kc = ',K_C)
    plt.show()