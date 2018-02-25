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
    """ 
    Main loop for calculating the membrene moduli, set the parameters and 
    gets the miduli for each frame. Stores the means result 
    """
    
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
    t_long=[] #lonitual director componants spectra
    t_trans=[] #tansverse director component spectra
    n_long_com=[]
    n_trans_com =[]
    s_wavenum = [] #flattned wavenumbers
    s_hf = [] # flatened height spectra
    
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
# =============================================================================
# Collecting the Director spectra
# =============================================================================
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
            
            n_qx=gf.fourier_trans_grid(director_grid_av[:,:,0],
                                       M,test_frame.box_size[0])
            n_qy=gf.fourier_trans_grid(director_grid_av[:,:,1],
                                       M,test_frame.box_size[0])
            n_q = np.zeros((n_qx[0,:].size, 
                            n_qx[:,0].size, 2), dtype=np.complex_)
            n_q[:,:,0] = n_qx
            n_q[:,:,1] = n_qy
            #collecting the spectra
            n_pow_sum, w_grid_sum, n_comp = gf.collect(n_q,q_grid,M,box_size[0])
            sumcoll = n_pow_sum #/frames_per_block
            gridcoll = w_grid_sum #/frames_per_block
            
# =============================================================================
# Collecting the height spectra
# =============================================================================
            height_field_up = gf.create_height_grid(lipid_grid_upper,M, 
                                               test_frame.box_size,
                                               test_frame.bilayer.mean_z)
            height_field_down = gf.create_height_grid(lipid_grid_lower,M, 
                                               test_frame.box_size,
                                               test_frame.bilayer.mean_z)
            hf_av = 0.5*(height_field_up  + height_field_down)
            hf_q = gf.fourier_trans_grid(hf_av,
                                       M,test_frame.box_size[0])
            hf_q *= (box_size[0]/(M**2))/10
            hf_pow = np.square(np.absolute(hf_q))
            
           
 # =============================================================================
# Collecting the tilt spectra
# =============================================================================         
            #nx_up,ny_up = gf.create_normal_grid(
            #        height_field_up, M,test_frame.box_size[0])
            #nx_down,ny_down = gf.create_normal_grid(
            #        height_field_down, M,test_frame.box_size[0])
            #nx_up,ny_up = np.gradient(height_field_up,test_frame.box_size[0]/M)
            #nx_down,ny_down = np.gradient(height_field_down,test_frame.box_size[0]/M)
            
            #nx_up *= M**2
            #ny_up *= M**2
            #nx_down *= M**2
            #ny_down *= M**2
            
            nx,ny = gf.create_normal_grid(
                    hf_av, M,test_frame.box_size[0])
            #nx, ny = np.gradient(hf_av,test_frame.box_size[0]/M)
# =============================================================================
#             tiltx_up = director_grid_upper[:,:,0] - nx_up
#             tilty_up = director_grid_upper[:,:,1] - ny_up
#             tiltx_down = director_grid_lower[:,:,0] - nx_down
#             tilty_down = director_grid_lower[:,:,1] - ny_down
# =============================================================================
            
            tiltx = nx #0.5*(nx_up + nx_down)
            tilty = ny #0.5*(ny_up + ny_down)
            tqx =nx # gf.fourier_trans_grid(tiltx,
#                                       M,test_frame.box_size[0])
            tqy=ny # gf.fourier_trans_grid(tilty,
    #                                  M,test_frame.box_size[0])
            tq = np.zeros((tqx[0,:].size,tqx[:,0].size,2), dtype=np.complex_)
            tq[:,:,0] = tqx
            tq[:,:,1] = tqy
            t_pow_sum, w_grid_sum2, t_comp = gf.collect(tq,q_grid,M,box_size[0])
            # and a 1D array containing the q values and correspondant logitudenal/
            # Transverse componants of the director field
            sorted_wavenum_u, n_long_u, n_trans_u, n_long_c, n_trans_c, hf_u, tl , tp = (
                    gf.get_moduli(gridcoll,sumcoll[:,:,0],sumcoll[:,:,1],
                                  n_comp[:,:,0], n_comp[:,:,1],hf_pow,
                                  t_pow_sum[:,:,0], t_pow_sum[:,:,1]))
            tl = sorted_wavenum_u**2*hf_u - n_trans_u
            n_long.append(n_long_u)
            n_trans.append(n_trans_u)
            t_long.append(tl)
            t_trans.append(tp)
            n_long_com.append(n_long_c)
            n_trans_com.append(n_trans_c)
            s_wavenum.append(sorted_wavenum_u)
            s_hf.append(hf_u)
            

            
            print ("\r" + "{0:.5f}".format(100*(i/frames_per_block)) 
                   + "% done \r",end="")      
    
    #converting to numpy matrices 
    trans = np.array(n_trans)
    long = np.array(n_long)
    ttrans = np.array(t_trans)
    tlong = np.array(t_long)
    long_c = np.array(n_long_com)
    trans_c = np.array(n_trans_com)
    hf = np.array(s_hf)
    
    save_dict = {'q': s_wavenum[0],
                 'n_long': long.mean(axis=0),
                 'n_long_se': long.std(axis=0)/np.sqrt(number_of_frames),
                 'n_trans': trans.mean(axis=0),
                 'n_tans_se': trans.std(axis=0)/np.sqrt(number_of_frames),
                 't_long': tlong.mean(axis=0),
                 't_long_se': tlong.std(axis=0)/np.sqrt(number_of_frames),
                 't_trans': ttrans.mean(axis=0),
                 't_tans_se': ttrans.std(axis=0)/np.sqrt(number_of_frames),
                 'hf':np.mean(hf, axis=0)}
    scipy.io.savemat(params.output+'.mat', mdict=save_dict)
    #creating and saving pandas dataframe
    
    df = pd.DataFrame(save_dict)
    df.set_index('q', inplace=True)
    writer = pd.ExcelWriter(params.output+'.xlsx')
    df.to_excel(writer, sheet_name='Director spectra')
    writer.save()
    df.to_csv(params.output+'.csv')
    
    # Saving numpy files
    np.save('trans{}'.format(params.lipids),trans)
    np.save('long{}'.format(params.lipids), long)
    np.save('trans_c{}'.format(params.lipids),trans_c)
    np.save('long_c{}'.format(params.lipids), long_c)
    
    np.save('wavenum', s_wavenum[0])
     
    # Just for testing, remove before final submission
 #=============================================================================
#    K_C_all = np.zeros(number_of_blocks)
#    for i in range(number_of_blocks):
#        K_C_tmp = np.square(s_wavenum)*trans
#        K_C_all[i] = k_b*T/np.mean(K_C_tmp[0:2])
#     
#    sorted_wavenum_u, n_long_u, n_trans_u, nlc, ntc = gf.get_moduli(
#                                            gridcoll,sumcoll[:,:,0],
#                                            sumcoll[:,:,1],
#                                            n_comp[:,:,0], n_comp[:,:,1],hf)
#    plt.plot(sorted_wavenum_u,n_long_u)
#    K_C = np.square(sorted_wavenum_u)*n_trans_u
#    K_C = k_b*T/np.mean(K_C[0:2])
#    print('Kc = ',K_C)
#    plt.show()
# =============================================================================
