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

fig308,ax308=plt.subplots(1,2)
fig108,ax108=plt.subplots(1,2)
fig305,ax305=plt.subplots(1,2)

ax308[0].plot(bright308[0],bright308[1],'x')
ax308[1].plot(bright308[2],bright308[3][1100,:],'-x')
ax308[0].set_title('1100308004 Brightness')
ax308[1].set_title('1100308004 Emissivity')

ax108[0].plot(bright108[0],bright108[1],'x')
ax108[1].plot(bright108[2],bright108[3][1100,:],'-x')
ax108[0].set_title('1080416025 Brightness')
ax108[1].set_title('1080416025 Emissivity')

ax305[0].plot([i+1 for i in range(len(bright305[1]))],bright305[1],'-x')
ax305[1].plot(bright305[2],bright305[3][1100,:],'-x')
ax305[0].set_title('1100305023 Brightness')
ax305[1].set_title('1100305023 Emissivity')

