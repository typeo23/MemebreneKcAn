#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 10:55:19 2017

@author: itay
Utilitie functions for creating a grid from lipid list and interpulating
"""
import numpy as np
from scipy import interpolate
from scipy.interpolate import griddata


def create_lipid_grid(lipids, M, L):
    """ Create a MxM grid from a lipid list lipids, box size is L """
    grid = [[[] for i in range(M)]for i in range(M)]
    dLx = L[0]/float(M)
    dLy = L[1]/float(M)
    for lipid in lipids:
        gridX = lipid.head[0]/float(dLx)
        gridY = lipid.head[1]/float(dLy)
        k_ind = int(np.floor(gridX))
        l_ind = int(np.floor(gridY))
        " Preiodic boudary conditions"
        if (k_ind < 0):
            k_ind += M
        if (k_ind >=M):
            k_ind -= M
        if (l_ind <0):
            l_ind += M
        if (l_ind >= M):
            l_ind -= M
        #direc = lipid.head - lipid.tail
        grid[k_ind][l_ind].append(lipid)
    return grid


                
    
# =============================================================================
# ===Commented for now, will be removed if we decide not to use closest point==
#    
# def get_closest_point(lipds,location,L):
#     """ return the closest lipid in lipds list to location"""
#     min_dist = L
#     for lipid in lipds:
#         lipid_xy = np.array([lipid.tail[0], lipid.tail[1]])
#         distance_from_gridpoint = np.linalg.norm(location - lipid_xy)
#         if (distance_from_gridpoint < min_dist):
#             min_dist = distance_from_gridpoint
#             closest_lipid = lipid;
#     return closest_lipid
# 
# 
# def create_lipid_grid_closest2(lipids,M,L,n,pool):    
#     """ Create MxM lipid grid using the closest lipid to the gridpoint
#     Accelareted by grdding the lipids and only scanning n surronding blocks"""
#     
#     lipid_grid = create_lipid_grid(lipids, M, L)
#     out_grid = [[[] for i in range(M)] for i in range(M)]
#     for k_ind in range(M):
#         for l_ind in range(M):
#             lipid_list = []
#             # creates a list of the lipids in the n surrounding blocks
#             for ind_x in range(k_ind-n, k_ind+n):
#                 for ind_y in range(l_ind-n, l_ind+n):
#                     currx = ind_x
#                     curry = ind_y
#                 
#                     if (currx <0):
#                         currx += M
#                     if (currx > M-1):
#                         currx -= M
#                     
#                     if (curry <0):
#                         curry += M
#                     if (curry > M-1):
#                         curry -= M
#                     for lipid in lipid_grid[currx][curry]:
#                         lipid_list.append(lipid)
#             # now get the colsest one
#             location = np.array([k_ind*L[0]/M, l_ind*L[0]/M])        
#             out_grid[k_ind][l_ind] = [get_closest_point(lipid_list,
#                                                         location,L[0])]
#     return out_grid
#             
#     
# =============================================================================
    
    
    
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
    return np.mean([lipid.dir for lipid in lipids], axis=0)

def create_height_grid(lipid_grid, M, L, z_mean):
    """ 
    Calculets the height field from the lipid_grid, z_mean is the mean
    height if the lipid heads. 
    //TODO: if we decide to keep it condier refactoring to all grid avraging 
    functions """
    #first calculate the hight field
    empty_grid_points =[]
    height_grid = np.zeros((M,M))
    for k_ind in range(M):
        for l_ind in range(M):
            if (len(lipid_grid[k_ind][l_ind])>0):
                height_grid[k_ind,l_ind] = np.mean(
                                        [lipid.head[2] for lipid in 
                                        lipid_grid[k_ind][l_ind]])- z_mean
            else:
                empty_grid_points.append((k_ind, l_ind))
                
    height_grid = interpulate_grid(height_grid, empty_grid_points, M)
    return height_grid
                
                
def create_normal_grid(height_grid,M,L):
    """
    Create the normal grid from the height grid (N = -iqh)
    """
    spacing = L / float(M)

    q1 = 2*np.pi*np.fft.fftfreq(M, d=spacing)
    q2 = 2*np.pi*np.fft.fftfreq(M, d=spacing)
    #q1 = 1j*np.zeros(M)
    #q1[0:int(M/2.0)] = np.arange(0, int(M/2.0), 1)
    #q1[int(M/2)+1:] = np.arange(int(-M/2.0) + 1, 0, 1)
    #q1 = (2*np.pi/L)*q1
    
    qx, qy = np.meshgrid(q2,q1)
#    hq = np.fft.fft2(height_grid)
#    norm_x = np.fft.ifft2(qx*hq)
#    norm_y = np.fft.ifft2(qy*hq)
    #normx = np.zeros((M,M), dtype=np.complex)
    #normy = np.zeros((M,M), dtype=np.complex)
    norm = 1 # L/M**2
    twoPiLx = (2*np.pi)/L
    hq = np.fft.fft2(height_grid)
    #hq2 = np.fft.np.fft.rfft2(np.transpose(height_grid))
    
# =============================================================================
#     for i in range(M):
#         for j in range(M):
#             if (i < M/2):
#                 qx = i
#             else:
#                 qx = i-M
#                 
#             if (j < M/2):
#                 qy = i
#             else:
#                 qy = i-M
#                 
#             if i == M/2:
#                 qx =0
#             if j == M/2:
#                 qy = 0
#             realx = -norm*q1[i]*hq[i][j].imag*twoPiLx
#             imagx = norm*q1[i]*hq[i][j].real*twoPiLx
#             normx[i][j] = realx + 1j*imagx
#             
#             realy = -norm*q1[j]*hq[i][j].imag*twoPiLx
#             imagy = norm*q1[j]*hq[i][j].real*twoPiLx
#             normy[i][j] = realy + 1j*imagy
# =============================================================================
            
    normx = 1j*qx*hq
    normy = 1j*qy*hq
    #norm_x = np.fft.irfft2(normx)
    #norm_y = np.transpose(np.fft.irfft2(normy))
           
    
    return normx, normy
    

def interpulate_grid(grid, empty_grid_points,M):
    """ bi-linear interpulation for emty grid points. Assumes there are no
    two adjucnt empty patches"""
    for empty_point in empty_grid_points:
        for ind_x in range(empty_point[0]-1, empty_point[0]+1):
            for ind_y in range(empty_point[1]-1, empty_point[1]+1):
                currx = ind_x
                curry = ind_y
                # periodic boundary conditions
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

def interpulate_grid2 (grid,empty_grid_points,M):
    """ Interpulates the grid with arbitrary missing points"""
    non_zero_indices =np.nonzero(grid[:,:,2])
    non_zero_vals = grid[np.nonzero(grid[:,:,2])]
    non_zero_vals = np.reshape(non_zero_vals ,
                               (int(non_zero_indices[0].size),3))
    x_grid_ind = non_zero_indices[0]
    y_grid_ind = non_zero_indices[1]
    
    x_val = non_zero_vals[:,0]
    y_val = non_zero_vals[:,1]
    z_val = non_zero_vals[:,2]
    
    grid_x, grid_y = np.mgrid[0:M, 0:M]
    
    grid_x_int = griddata(non_zero_indices, x_val, 
                          (grid_x,grid_y), fill_value=0,method='linear')
    
    grid_y_int = griddata(non_zero_indices, y_val, 
                          (grid_x,grid_y), fill_value=0,method='linear')
    
    grid_z_int = griddata(non_zero_indices, z_val, 
                          (grid_x,grid_y), fill_value=0,method='linear')
    grid_int = np.dstack((grid_x_int,grid_y_int,grid_z_int))
   
    grid2=interpulate_grid(grid, empty_grid_points,M) 
    bb = np.multiply(grid==0,grid2)
    return (grid_int)

def normalize_grid(grid,M):
    """  "Nomalize" the vector values at each gridpoint, such that Z=1"""
    for ind_x in range(M):
        for ind_y in range(M):
            if grid[ind_x, ind_y, 2] != 0:
                grid[ind_x, ind_y] /= grid[ind_x, ind_y, 2]
               # grid[ind_x, ind_y] /= np.sqrt(
                #        grid[ind_x, ind_y, 0]**2+grid[ind_x, ind_y, 1]**2
                #        + grid[ind_x, ind_y, 2]**2)
                

def create_grid_qvalues(M, L):
    """ return a MxM grid (of size L in real space) 
    of fourier wave vectors q"""  
    
    spacing = L[0] / float(M)

    wave_num1 = 2*np.pi*np.fft.fftfreq(M, d=spacing)
    wave_num2 = 2*np.pi*np.fft.fftfreq(M, d=spacing)
    
# =============================================================================
#     q1 = np.zeros(M)
#     q1[0:int(M/2.0)] = np.arange(0, int(M/2.0), 1)
#     q1[int(M/2)+1:] = np.arange(int(-M/2.0) + 1, 0, 1)
#     wave_num1 = (2*np.pi/L[0])*q1
#     wave_num2 = (2*np.pi/L[0])*q1
# =============================================================================
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
    Remove the Nyquist components as they diverge
    """
    n_q1 = np.fft.fft2(grid_val_real, norm = None )
    
    # Remove areas of the grids, which correspond to wavenumbers with 
    #components corresponding to the Nyquist frequencies    
    if ( M % 2 == 0):
        n_q1 = np.delete(n_q1, M / 2, axis = 0)
        n_q1 = np.delete(n_q1, M / 2, axis = 1)
 

    else:   # This statement should not be necessary as M should be even
        n_q1 = np.delete(n_q1, (M + 1) / 2, axis = 0)
        n_q1 = np.delete(n_q1, (M + 1) / 2, axis = 1)
        n_q1 = np.delete(n_q1, (M - 1) / 2, axis = 0)
        n_q1 = np.delete(n_q1, (M - 1) / 2, axis = 1)
                                         
    return n_q1

def collect(n_q, w_grid, M,L):
    """
    Function to calculate the orientation fields of each monolayer, 
    then perform 2D FFTs on the data sets. Power spectra are then calculated 
    from the Fourier data and stored in a variable containing the sum of 
    the spectra over all sampled steps
    """
    
    n_q *= (L/(M**2))   # Scaling factor to correct for units and account for grid size dependency

    # DFT for  the field, and seperate the longitodenal 
    #and transverse components
    siz1 = w_grid[:, 0, 0].size
    siz2 = w_grid[0, :, 0].size
    q_matrix = np.zeros((2, 2))
    n_comp = np.zeros((siz1, siz2, 2), dtype = np.complex_)
    for i in range(siz1):
        for j in range(siz2):
            if ( i !=0 or j != 0):
                q_v = w_grid[i, j, :]
                n_v = n_q[i, j, :]

                #q = np.linalg.norm( q_v ) # Magnitude of wavenumber
                q = np.sqrt(q_v[0]**2 + q_v[1]**2)
                q_matrix[0,:] = q_v
                q_matrix[1,0] = - q_v[1]
                q_matrix[1,1] = q_v[0]

                n_comp[i,j,:] = np.dot(q_matrix, n_v) / q

    # Calculate power spectra of n_trans and n_long
    n_pow = np.square(np.absolute(n_comp)) # Contains power spectrum of n_trans and n_long
    n_pow_sum =  n_pow
    w_grid_sum =  w_grid
    return n_pow_sum, w_grid_sum, n_comp


def get_moduli(w_grid_av, n_trans_av, 
               n_long_av, n_trans_c, n_long_c, height_f, t_long_av, t_trans_av):
    """
    return a flattened set of unique wavenumbers and the corresponding grid 
    values for the longitudenal and transverse components
    """
    w_grid_mag = np.linalg.norm(w_grid_av, axis = 2)

    # First flatten the arrays and remove the q=(0,0) wavenumber values
    n_trans_av = n_trans_av.flatten()
    n_trans_av= n_trans_av[1:]
    n_long_av = n_long_av.flatten()
    n_long_av = n_long_av[1:]
    t_trans_av = t_trans_av.flatten()
    t_trans_av= t_trans_av[1:]
    t_long_av = t_long_av.flatten()
    t_long_av = t_long_av[1:]
    n_trans_c = n_trans_c.flatten()[1:]
    n_long_c = n_long_c.flatten()[1:]
    w_grid_mag = w_grid_mag.flatten()
    w_grid_mag=w_grid_mag[1:]
    height_f = height_f.flatten()
    height_f = height_f[1:]
	
    # Sort the wavenumbers and then accordingly rearrange the other arrays
    sorted_indices = np.argsort(w_grid_mag)
    sorted_wavnum = w_grid_mag[sorted_indices]
    sorted_n_trans = n_trans_av[sorted_indices]
    sorted_n_long = n_long_av[sorted_indices]
    sorted_t_trans = t_trans_av[sorted_indices]
    sorted_t_long = t_long_av[sorted_indices]
    sorted_n_trans_c = n_trans_c[sorted_indices]
    sorted_n_long_c = n_long_c[sorted_indices]
    sorted_heigh_f = height_f[sorted_indices]

    # Take the average of the data for each unique wavenumber

    # Returns the unique wavenumber and the indices of the first 
    #element containing each unique wavenumber
    sorted_wavnum_u, unique_ind = np.unique(sorted_wavnum, return_index=True)
    cnt = 0
    n_trans_u = np.zeros(sorted_wavnum_u.size)
    n_long_u = np.zeros(sorted_wavnum_u.size)
    t_trans_u = np.zeros(sorted_wavnum_u.size)
    t_long_u = np.zeros(sorted_wavnum_u.size)
    n_trans_c = np.zeros(sorted_wavnum_u.size, dtype = np.complex_)
    n_tlong_c = np.zeros(sorted_wavnum_u.size, dtype = np.complex_)
    height_f_u = np.zeros(sorted_wavnum_u.size)
    for index in unique_ind:
        # Get boolean array of elements in which true elements correspond 
        #to the unique wavenumber
        condlist = (sorted_wavnum == sorted_wavnum[index])

        # Get all n_trans/n_long with a given unique wavenumber and 
        #calc their average values
        n_trans = np.compress(condlist, sorted_n_trans)
        n_trans_u[cnt] = np.mean(n_trans)

        n_long = np.compress(condlist, sorted_n_long)
        n_long_u[cnt] = np.mean(n_long)
        
        t_trans = np.compress(condlist, sorted_t_trans)
        t_trans_u[cnt] = np.mean(t_trans)

        t_long = np.compress(condlist, sorted_t_long)
        t_long_u[cnt] = np.mean(t_long)
        
        n_transc = np.compress(condlist, sorted_n_trans_c)
        n_trans_c[cnt] = np.mean(n_transc)

        n_longc = np.compress(condlist, sorted_n_long_c)
        n_long_c[cnt] = np.mean(n_longc)
        
        height_f = np.compress(condlist,  sorted_heigh_f)
        height_f_u[cnt] = np.mean(height_f)

        cnt +=1

    # Convert units from Angstroms to nm
    sorted_wavnum_u *= 10  # Angstrom^-1 -> nm^-1
    n_long_u *= 0.01       # Angstrom^2  -> nm^2
    n_trans_u *= 0.01      # Angstrom^2  -> nm^2
    t_long_u *= 0.01       # Angstrom^2  -> nm^2
    t_trans_u *= 0.01      # Angstrom^2  -> nm^2
    n_long_c *= 0.01
    n_trans_c *= 0.01
    height_f_u *= 0.01

    return (sorted_wavnum_u, n_long_u, 
            n_trans_u, n_trans_c, n_long_c, height_f_u, t_long_u, t_trans_u)


