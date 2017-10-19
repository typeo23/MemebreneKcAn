#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 18:34:37 2017

@author: itay
A full simulation frame, contains a lipid bilayer and box dimentions
Also contains routines to read a frame from a file
"""
import numpy as np
from Lipid_bilayer  import Lipid_bilayer as lipids
from Lipid_bilayer  import Lipid as lipid

class Frame():
    def __init__(self,lipids_number,lipidx='lipidX.out',
                 lipidy='lipidY.out',lipidz='lipidZ.out'
                 ,boxx='boxSizeX.out',boxy='boxSizeY.out',boxz='boxSizeZ.out'):
        """ Initialize a frame containing a lipid bilayer and box dimentions 
        will be loaded from the files lipidx,y,z and box dimentiond will be 
        loaded from boxx,y,z """
        self.lipids_number = lipids_number
        self.lipids_x_file = open(lipidx,'r')
        self.lipids_y_file = open(lipidy,'r')
        self.lipids_z_file = open(lipidz,'r')
        self.box_x_file = open(boxx,'r')
        self.box_y_file = open(boxy,'r')
        self.box_z_file = open(boxz,'r')
        
    def load_next_frame(self):
        """ read the next frame of lipi bilayer """
        bilayer = [];
        self.box_size =[]
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
            
#            for i in range(3):
#                if (head[i] < 0):
#                    head[i] += self.box_size[i]
#                if (head[i] >= self.box_size[i]):
#                    head[i] -= self.box_size[i]
           
            tail1.append(float(self.lipids_x_file.readline()))
            tail1.append(float(self.lipids_y_file.readline())) 
            tail1.append(float(self.lipids_z_file.readline())) 
#            for i in range(3):
#                if (tail1[i] < 0):
#                    tail1[i] += self.box_size[i]
#                if (tail1[i] >= self.box_size[i]):
#                    tail1[i] -= self.box_size[i]
                    
            tail=np.array(tail1)
            
#            tail2.append(float(self.lipids_x_file.readline()))
#            tail2.append(float(self.lipids_y_file.readline())) 
#            tail2.append(float(self.lipids_z_file.readline())) 
#            
#            for i in range(3):
#                if (tail2[i] < 0):
#                    tail2[i] += self.box_size[i]
#                if (tail2[i] >= self.box_size[i]):
#                    tail2[i] -= self.box_size[i]
#            
#            tail2=np.array(tail2)
#            
#            tail = 0.5*(tail1+tail2)
            
            bilayer.append(lipid(head,tail))
        
        self.bilayer = lipids(bilayer)
            