# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:59:41 2019

@author: rmreksoatmodjo

Collection of general Tools to perform oft-repeated SOLPS data analyis and post-processing tasks
"""
import os
import numpy as np

def SET_WDIR(BASEDRT,TOPDRT): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'rmreksoatmodjo':
            BASEDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/rmreksoatmodjo/Desktop/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == '18313':
            BASEDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"G:/My Drive/College of William and Mary/Research/SOLPS Stuff/"
        elif os.environ['USERNAME'] == 'Richard':
            BASEDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/SOLPS_2D_prof/"
            TOPDRT = r"C:/Users/Richard/WMGDrive/College of William and Mary/Research/SOLPS Stuff/"
    else:
        BASEDRT=BASEDRT
        TOPDRT=TOPDRT
    
    return BASEDRT, TOPDRT

def TANH(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def EXPFIT(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)
'''
class FORT44(object): #Class of methods used to parse and organize SOLPS fort.44 data
    
    def __init__(self,Shot,Attempt,Parameter,**kwargs):
        
        self._reset_object()
        
        self.Shot=Shot
        
        self.Attempt=Attempt
        
        self.Parameter=Parameter
        
    def _reset_object(self):
        self.Shot=None
        self.Attempts=None
        self.Parameter=None
        
    def _OpenFort44(self):
        
    def _FindParameter(self):
        
    def PlotWalls(self):
'''        
        