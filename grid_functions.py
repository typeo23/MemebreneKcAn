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
            
    
