# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:05:43 2019

@author: rmreksoatmodjo
"""

import xarray as xr
from VesselPlotter import SOLPSPlotter
import pandas as pd
import numpy as np

JXA = 55
Shot = 25
Attempt = 129
PARAM = 'DN'
Title = 'Shot025_D'

ShotData = SOLPSPlotter(Shot,[Attempt],PARAM,'Export',Publish=[],SEP=20,LOG10=0)

RR = ShotData['RR'].loc[:,JXA,Attempt].to_dataframe()
PRM = ShotData['PARAM'].loc[:,JXA,Attempt].to_dataframe()

PRM_RR = pd.concat([RR,PRM],axis=1)

PRM_RR.to_csv(Title,columns=(RR.columns[2],PRM.columns[2]))

  