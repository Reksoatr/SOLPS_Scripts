# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:59:41 2019

@author: rmreksoatmodjo

Function to set correct Working Directory Path depending on which machine is in use
"""
import os

def SET_WDIR(BASEDRT,TOPDRT):
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'rmreksoatmodjo':
            BASEDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == '18313':
            BASEDRT = r"C:/Users/18313/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/18313/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Richard':
            BASEDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
    else:
        BASEDRT=BASEDRT
        TOPDRT=TOPDRT
    
    return BASEDRT, TOPDRT