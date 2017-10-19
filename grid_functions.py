#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 10:55:19 2017

@author: itay
Utilitie functions for creating a grid from lipid list and interpulating
"""
import numpy as np


def create_lipid_grid(lipids, M, L):
    """ Create a MxM grid from a lipid list lipids, box size is L """
    grid = [[[] for i in range(M)]for i in range(M)]
    dLx = L[0]/float(M)
    dLy = L[1]/float(M)
    for lipid in lipids:
        gridX = lipid.tail[0]/float(dLx)
        gridY = lipid.tail[1]/float(dLy)
        k_ind = int(np.floor(gridX))
        l_ind = int(np.floor(gridY))
        if (k_ind < 0):
            k_ind += M
        if (k_ind >=M):
            k_ind -= M
        if (l_ind <0):
            l_ind += M
        if (l_ind >= M):
            l_ind -= M
        grid[k_ind][l_ind].append(lipid)
    return grid


def create_director_grid(lipid_grid, M, L):
    """ Create a numpy 2d MxM grid of the avrage direcotr from a lipid grid
    Interpulate empty patches
    """
    director_grid = np.zeros((M, M, 3))
    empty_grid_points =[] 
    
    for k_ind in range(M):
        for l_ind in range(M):
            if (len(lipid_grid[k_ind][l_ind]) > 0):
                director_grid[k_ind][l_ind][:] = calc_director(lipid_grid[k_ind][l_ind])
            else:
                empty_grid_points.append((k_ind, l_ind))
    return interpulate_grid(director_grid,empty_grid_points,M)

def calc_director(lipids):
    """ calculaet the avrage director from a list of lipids in a patch """
    lipid_num = len(lipids)
    av_director = np.zeros(3)
    for lipid in lipids:
        av_director += lipid.dir
    av_director /= lipid_num
    return av_director

def interpulate_grid(grid, empty_grid_points,M):
    for empty_point in empty_grid_points:
        for ind_x in range(empty_point[0]-1, empty_point[0]+1):
            for ind_y in range(empty_point[1]-1, empty_point[1]+1):
                currx = ind_x
                curry = ind_y
                
                if (currx <0):
                    currx += M
                if (currx > M-1):
                    currx -= M
                    
                if (curry <0):
                    curry += M
                if (curry > M-1):
                    curry -= M
                grid[empty_point] += grid[currx, curry]
                
        grid[empty_point] /= 8
    return grid


def normalize_grid(grid,M):
    for ind_x in range(M):
        for ind_y in range(M):
            grid[ind_x, ind_y] /= grid[ind_x, ind_y, 2]
   

def create_grid_qvalues(M, L):
    """ return a MxM grid (of size L in rela space) 
    of fourier wave vectors q"""  
    
    spacing = L[0] / float(M)

    wave_num1 = 2*np.pi*np.fft.fftfreq(M, d=spacing)
    wave_num2 = 2*np.pi*np.fft.fftfreq(M, d=spacing)

    # Remove the Nyquist frequencies
    if (M % 2 == 0):
        wave_num1 = np.delete(wave_num1, M / 2, None)
        wave_num2 = np.delete(wave_num2, M / 2, None)
    else:  # Should not be necessary as  should be even
        wave_num1 = np.delete(wave_num1, (M + 1) / 2, None)
        wave_num2 = np.delete(wave_num2, (M + 1) / 2, None)
        wave_num1 = np.delete(wave_num1, (M - 1) / 2, None)
        wave_num2 = np.delete(wave_num2, (M - 1) / 2, None)

    # Create an MxM grid of 2D wavenumbers
    siz1 = wave_num1.size
    siz2 = wave_num2.size

    w_grid = np.zeros((siz1, siz2, 2))
    for i in range(siz1):
        for j in range(siz2):
            w_grid[i, j, 0] = wave_num1[i]
            w_grid[i, j, 1] = wave_num2[j]

    return w_grid

def fourier_trans_grid(grid_val_real, M, L):
    """
    Calculate the 2D FFT of grid_val_real.
    FFTs are carried out separately on the x- and y-components of grid_val_real. 
    """
    n_q1 = np.fft.fft2( grid_val_real[:,:,0], norm = None )
    n_q2 = np.fft.fft2( grid_val_real[:,:,1], norm = None )

    # Remove areas of the grids, which correspond to wavenumbers with components corresponding to the Nyquist frequencies    
    if ( M % 2 == 0):
        n_q1 = np.delete(n_q1, M / 2, axis = 0)
        n_q1 = np.delete(n_q1, M / 2, axis = 1)

        n_q2 = np.delete(n_q2, M / 2, axis = 0)
        n_q2 = np.delete(n_q2, M / 2, axis = 1)  

    else:   # This statement should not be necessary as M should be even
        n_q1 = np.delete(n_q1, (M + 1) / 2, axis = 0)
        n_q1 = np.delete(n_q1, (M + 1) / 2, axis = 1)
        n_q1 = np.delete(n_q1, (M - 1) / 2, axis = 0)
        n_q1 = np.delete(n_q1, (M - 1) / 2, axis = 1)

        n_q2 = np.delete(n_q2, (M + 1) / 2, axis = 0)
        n_q2 = np.delete(n_q2, (M + 1) / 2, axis = 1)
        n_q2 = np.delete(n_q2, (M - 1) / 2, axis = 0)
        n_q2 = np.delete(n_q2, (M - 1) / 2, axis = 1)

    siz1 = len(n_q1[0,:])
    siz2 = len(n_q1[:,0])

    n_q = np.zeros((siz1, siz2, 2), dtype=np.complex_)
    n_q[:,:,0] = n_q1
    n_q[:,:,1] = n_q2
                                         
    return n_q

def collect(n_q, w_grid, M,L):
    """
    Function to calculate the orientation fields of each monolayer, 
    then perform 2D FFTs on the data sets. Power spectra are then calculated 
    from the Fourier data and stored in a variable containing the sum of 
    the spectra over all sampled steps
    """
    
    n_q *= (L/(M**2))   # Scaling factor to correct for units and account for grid size dependency

    # Transform n_q according to the SI of the 2014 paper by Levine et al. to get the
    # longitudinal and transverse components of the bilayer orientation fields
    siz1 = w_grid[:, 0, 0].size
    siz2 = w_grid[0, :, 0].size
    q_matrix = np.zeros((2, 2))
    n_comp = np.zeros((siz1, siz2, 2), dtype = np.complex_)
    for i in range(siz1):
        for j in range(siz2):
            if ( i !=0 or j != 0):
                q_v = w_grid[i, j, :]
                n_v = n_q[i, j, :]

                q = np.linalg.norm( q_v ) # Magnitude of wavenumber

                q_matrix[0,:] = q_v
                q_matrix[1,0] = - q_v[1]
                q_matrix[1,1] = q_v[0]

                n_comp[i,j,:] = np.dot(q_matrix, n_v) / q

    # Calculate power spectra of n_trans and n_long
    n_pow = np.square(np.absolute(n_comp)) # Contains power spectrum of n_trans and n_long
    n_pow_sum =  n_pow
    w_grid_sum =  w_grid
    return n_pow_sum, w_grid_sum

def get_moduli( w_grid_av, n_trans_av, n_long_av, plotg = False ):
    """
    The moduli are determined according to the model by Levine et al., 
    based on spectral analysis of lipid orientations.
    """
    k_b = 1.38064852e-23 # Boltzmann constant	
    w_grid_mag = np.linalg.norm(w_grid_av, axis = 2)

    # First flatten the arrays and remove the q=(0,0) wavenumber values
    n_trans_av = n_trans_av.flatten()
    n_trans_av= n_trans_av[1:]
    n_long_av = n_long_av.flatten()
    n_long_av = n_long_av[1:]
    w_grid_mag = w_grid_mag.flatten()
    w_grid_mag=w_grid_mag[1:]
	
    # Sort the wavenumbers and then accordingly rearrange the other arrays
    sorted_indices = np.argsort(w_grid_mag)
    sorted_wavnum = w_grid_mag[sorted_indices]
    sorted_n_trans = n_trans_av[sorted_indices]
    sorted_n_long = n_long_av[sorted_indices]

    # Take the average of the data for each unique wavenumber

    # Returns the unique wavenumber and the indices of the first element containing each unique wavenumber
    sorted_wavnum_u, unique_ind = np.unique(sorted_wavnum, return_index=True)
    cnt = 0
    n_trans_u = np.zeros(sorted_wavnum_u.size)
    n_long_u = np.zeros(sorted_wavnum_u.size)
    for index in unique_ind:
        # Get boolean array of elements in which true elements correspond to the unique wavenumber
        condlist = (sorted_wavnum == sorted_wavnum[index])

        # Get all n_trans/n_long with a given unique wavenumber and calc their average values
        n_trans = np.compress(condlist, sorted_n_trans)
        n_trans_u[cnt] = np.mean(n_trans)

        n_long = np.compress(condlist, sorted_n_long)
        n_long_u[cnt] = np.mean(n_long)

        cnt +=1

    # Convert units from Angstroms to nm
    sorted_wavnum_u *= 10  # Angstrom^-1 -> nm^-1
    n_long_u *= 0.01       # Angstrom^2  -> nm^2
    n_trans_u *= 0.01      # Angstrom^2  -> nm^2

    q_long_sq = np.square(sorted_wavnum_u)*n_long_u

    
    

    return sorted_wavnum_u, n_long_u, n_trans_u

