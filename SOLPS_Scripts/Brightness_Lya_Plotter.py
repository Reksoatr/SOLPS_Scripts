# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 00:17:59 2021

@author: 18313
"""

import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from TOOLS import SET_WDIR

BASE,TOP = SET_WDIR('','')

bright308=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lya_brightness_1100308004.pkl'.format(TOP),'rb'))
bright108=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lya_brightness_1080416025.pkl'.format(TOP),'rb'))
bright305=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lya_brightness_1100305023.pkl'.format(TOP),'rb'))

lyman308=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lyman_data_1100308004.pkl'.format(TOP),'rb'))
lyman108=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lyman_data_1080416025.pkl'.format(TOP),'rb'))
lyman305=pkl.load(open('{}gfileProcessing/cmod_files/Brightness_Profiles/lyman_data_1100305023.pkl'.format(TOP),'rb'))

plt.plot(bright308[0],bright308[1],'x')
plt.plot(bright308[2],bright308[3][1100,:],'x')