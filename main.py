#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Memebrene moduli calculations
Main program file

@author: itay
"""
import argparse 
from memmoduli import clac_moduli

def parse_command_line():
    """ 
    Add parameter to argparsers and return aggparse object 
    """
    
                          
    parse = argparse.ArgumentParser()
    # Input Files
    parse.add_argument("-lx", "--lipidx", 
                     help="A File with the lipids X positions", 
                     default="LipdX.out")
    parse.add_argument("-ly", "--lipidy", 
                     help="A File with the lipids Y positions", 
                     default="LipdY.out")
    parse.add_argument("-lz", "--lipidz", 
                     help="A File with the lipids X positions", 
                     default="LipdX.out")
    parse.add_argument("-bx", "--boxx", 
                     help="A File with the X axis box size", 
                     default="boxsizeX.out")
    parse.add_argument("-by", "--boxy", 
                     help="A File with the Y axis box size", 
                     default="boxsizeY.out")
    parse.add_argument("-bz", "--boxz", 
                     help="A File with the Z axis box size", 
                     default="boxsizeZ.out")
    #Simulation parameters
    parse.add_argument("-f", "--frames", 
                     help="Number of frmaes", type=int, required=True)
    parse.add_argument("-l", "--lipids", 
                     help="Number of lipids", type=int, required=True)
    parse.add_argument("-M", 
                     help="grid size", type=int, required=True)
    #output file
    parse.add_argument("-o", "--output", 
                     help="Outpu file name",  default="out")
     
    return parse
 
    
if (__name__ == "__main__"):
    """ parse argumants and calls the main program bodey"""
    args = parse_command_line()
    params= args.parse_args()
    clac_moduli(params)
      
 
 
