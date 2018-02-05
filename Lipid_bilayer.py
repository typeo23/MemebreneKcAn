# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 20:51:52 2017

@author: itay

Objects representing a single lipis and a lipid bilayer
"""

import numpy as np

class Lipid():
    def __init__(self,head,tail):
        self.head = head
        self.tail = tail
        self.dir = tail-head
    def __rpr__(self):
        return "<single lipid, contains head and tail vectors and director>"

class Lipid_bilayer():
    """ Class which contains all of the lipid information for a lipid bylayer 
    within a frame. """
    def __init__(self,lipids):
        """ initialize the bilayer with a list of lipids. The mean Z value
        is calculated and the lipids are assigned to up and down layers"""
        self.lipids = lipids
        self.mean_z = 0;
        self.upper = [];
        self.lower=[];
        self.calc_mean_z()
        self.devide_layers()
    def calc_mean_z(self):
        """ claculate mean Z value of the lipids head. Used to devide the lipid system 
        to upper and lower layer"""
        for lipid in self.lipids:
            self.mean_z += lipid.head[2]
        self.mean_z /= len(self.lipids)
    def devide_layers(self):
        for lipid in self.lipids:
            if(lipid.head[2] > self.mean_z):
                if(lipid.head[2]-lipid.tail[2] > 0):
                    self.upper.append(lipid)
            else:
                if (lipid.head[2]-lipid.tail[2] < 0):
                    self.lower.append(lipid)
    def __repr__(self):
        st = "<lipid bilayer containig %d lipids>" % len(self.lipids)
        return st 