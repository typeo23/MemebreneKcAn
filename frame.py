#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 18:34:37 2017

@author: itay
A full simulation frame, contains a lipid bilayer and box dimentions
Also contains routines to read a frame from a file
"""
import numpy as np
from Lipid_bilayer  import LipidBilayer as lipids
from Lipid_bilayer  import Lipid as lipid


class Frame():
    def __init__(self,lipids_number,lipidx='LipidX.out',
                 lipidy='LipidY.out',lipidz='LipidZ.out'
                 ,boxx='boxsizeX.out',boxy='boxsizeY.out',boxz='boxsizeZ.out'):
        """ 
        Initialize a frame containing a lipid bilayer and box dimentions 
        will be loaded from the files lipidx,y,z and box dimentiond will be 
        loaded from boxx,y,z 
        """
        self.lipids_number = lipids_number
        self.lipids_x_file = open(lipidx,'r')
        self.lipids_y_file = open(lipidy,'r')
        self.lipids_z_file = open(lipidz,'r')
        self.box_x_file = open(boxx,'r')
        self.box_y_file = open(boxy,'r')
        self.box_z_file = open(boxz,'r')
        self.box_size = []
        self.bilayer = None
        
    def load_next_frame(self):
        """ read the next frame of lipids from the file  """
        bilayer = []
        self.box_size.append(float(self.box_x_file.readline()))
        self.box_size.append(float(self.box_y_file.readline()))
        self.box_size.append(float(self.box_z_file.readline()))

        for lipid_num in range(self.lipids_number):
            head = []
            tail1 =[]
            tail2 = []
            head.append(float(self.lipids_x_file.readline())) 
            head.append(float(self.lipids_y_file.readline())) 
            head.append(float(self.lipids_z_file.readline())) 
            head = np.array(head)
            
#           
           
            tail1.append(float(self.lipids_x_file.readline()))
            tail1.append(float(self.lipids_y_file.readline())) 
            tail1.append(float(self.lipids_z_file.readline())) 
#            
                    
            tail = np.array(tail1)

            # filtering high angle lipids 
            dirc = head-tail
            dirc /= np.sqrt(dirc[0]**2 + dirc[1]**2 + dirc[2]**2)
            if np.fabs(dirc[2] > 0.5):
                bilayer.append(lipid(head, tail))
        
        self.bilayer = lipids(bilayer)
        
    def __del__(self):
        self.lipids_x_file.close()
        self.lipids_y_file.close()
        self.lipids_z_file.close()
        self.box_x_file.close()
        self.box_z_file.close()
