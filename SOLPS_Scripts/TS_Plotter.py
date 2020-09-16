# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:01:04 2020

@author: 18313
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider,TextBox,Button
from TOOLS import SET_WDIR
import numpy as np
import pickle as pkl
import sys

BASEDRT, TOPDRT=SET_WDIR('','')

Device='cmod'

Shots=['1120917011']#['1160718012','1160718013','1160718023']#['1160718024','1160718025']#['1160718012','1160718013','1160718023']#['1120917011']#['1101014006','1101014019']#    

Rad='psin'

Time0=75
AVG=0
PsinOffset=-0.005
TimeA=np.nan
TimeB=np.nan
CoreTS=True

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(3,1,height_ratios=[3,3,1],hspace=0.2)

ne_profile = fig.add_subplot(gs[0,0])
te_profile = fig.add_subplot(gs[1,0],sharex=ne_profile)
time_slide = fig.add_subplot(gs[2,0])

textTimeAax = plt.axes([0.15, 0.05, 0.03, 0.05])
TimeA_Text = TextBox(textTimeAax, 'Initial\nTime', hovercolor='0.9')

textTimeBax = plt.axes([0.2, 0.05, 0.03, 0.05])
TimeB_Text = TextBox(textTimeBax, 'End\nTime', hovercolor='0.9')

buttonTimeAVGax = plt.axes([0.25, 0.05, 0.11, 0.05])
TimeAVGButton = Button(buttonTimeAVGax, 'Plot Time Averaged Profiles')

N = len(Shots)
Data={}
NeLine={}
TeLine={}
Coords={'rmid':'(m)','psin':''}

for i in Shots:
    Data[i]=loadmat('{}gfileProcessing/{}_files/{}.mat'.format(TOPDRT,Device,i))
    if Rad not in Data[i].keys():
        sys.exit('Error! Coordinate {} does not exist! Aborting!'.format(Rad))
    if PsinOffset!=0:
        Data[i]['psin']=Data[i]['psin']+PsinOffset
    
    Data[i]['time']=Data[i]['time'].flatten()
    NeLine[i] = ne_profile.errorbar(Data[i][Rad][:,Time0],Data[i]['ne'][:,Time0],yerr=Data[i]['nerr'][:,Time0],marker='.',linestyle=':')
    TeLine[i] = te_profile.errorbar(Data[i][Rad][:,Time0],Data[i]['te'][:,Time0],yerr=Data[i]['terr'][:,Time0],marker='.',linestyle=':')
 
if CoreTS is True:
    with open('{}gfileProcessing/cmod_files/{}_CORE.pkl'.format(TOPDRT,Shots[0]),'rb') as f:
        CTS=pkl.load(f)
        
    CoreNe = ne_profile.errorbar(CTS[0],CTS[1]*1e20,yerr=CTS[2]*1e20,linestyle='',capsize=5,marker='.')
    CoreTe = te_profile.errorbar(CTS[3],CTS[4]*1000,yerr=CTS[5]*1000,linestyle='',capsize=5,marker='.')
    
ne_profile.set_title('Thompson Scattering Profiles at {:0.3f} sec'.format(Data[i]['time'][Time0]))
ne_profile.set_ylabel(r'Electron Density $n_e\;(m^{-3})$')
ne_profile.legend(Shots)
ne_profile.axhline(0.0,color='k')
ne_profile.grid()

te_profile.set_ylabel(r'Electron Temperature $T_e\;(eV)$')
te_profile.set_xlabel('{} {}'.format(Rad,Coords[Rad]))
te_profile.axhline(0.0,color='k')
te_profile.grid()

def update(event):
   global TimeA
   global TimeB
   global AVG

   for i in Shots:
       if AVG==1:
           print('Taking Median Average between {} and {} seconds'.format(TimeA,TimeB))
           ti = 0
           while float(TimeA) > Data[i]['time'].flatten()[ti]:
               ti = ti+1           
       
           tf = 0
           while float(TimeB) > Data[i]['time'].flatten()[tf]:
               tf = tf+1
           
           x_base_ne=np.median(Data[i][Rad][:,ti:tf],axis=1)
           y_base_ne=np.median(Data[i]['ne'][:,ti:tf],axis=1)
           y_error_ne=np.median(Data[i]['nerr'][:,ti:tf],axis=1)
           
           x_base_te=np.median(Data[i][Rad][:,ti:tf],axis=1)
           y_base_te=np.median(Data[i]['te'][:,ti:tf],axis=1)
           y_error_te=np.median(Data[i]['terr'][:,ti:tf],axis=1)
           
           TitleText='{:0.3f} to {:0.3f}'.format(Data[i]['time'].flatten()[ti],Data[i]['time'].flatten()[tf])
           
       else:
       
           Time1 = int(TimeSlider.val)
        
           x_base_ne=Data[i][Rad][:,Time1]
           y_base_ne=Data[i]['ne'][:,Time1]
           y_error_ne=Data[i]['nerr'][:,Time1]
           
           x_base_te=Data[i][Rad][:,Time1]
           y_base_te=Data[i]['te'][:,Time1]
           y_error_te=Data[i]['terr'][:,Time1]
           
           TitleText='{:0.3f}'.format(Data[i]['time'][Time1])
           
           TimeSlider.valtext.set_text('{:0.3f}\nseconds'.format(Data[i]['time'][Time1]))

       yerr_top_ne = y_base_ne + y_error_ne
       yerr_bot_ne = y_base_ne - y_error_ne
       
       new_segments_y_ne = [np.array([[x, yt], [x,yb]]) for x, yt, yb in zip(x_base_ne, yerr_top_ne, yerr_bot_ne)]
       
       NeLine[i][0].set_data(x_base_ne,y_base_ne)
       NeLine[i][2][0].set_segments(new_segments_y_ne)

       yerr_top_te = y_base_te + y_error_te
       yerr_bot_te = y_base_te - y_error_te
       
       new_segments_y_te = [np.array([[x, yt], [x,yb]]) for x, yt, yb in zip(x_base_te, yerr_top_te, yerr_bot_te)]
       
       TeLine[i][0].set_data(x_base_te,y_base_te)
       TeLine[i][2][0].set_segments(new_segments_y_te)
   
   ne_profile.relim()
   ne_profile.autoscale_view()
   ne_profile.set_title('Thompson Scattering Profiles at {} seconds'.format(TitleText))    
   te_profile.relim()
   te_profile.autoscale_view()
   
   fig.canvas.draw_idle()

TimeSlider = Slider(time_slide, 'Time', 0, len(Data[i]['time']), valinit=Time0, valfmt="%i", valstep=1)
TimeSlider.valtext.set_text('{:0.3f}\nseconds'.format(Data[i]['time'][Time0]))
TimeSlider.on_changed(update)

def submitTimeA(text):
    global TimeA
    if text == '':
        TimeA=np.nan
    else:
        TimeA=text    
    print('Initial Time={}'.format(TimeA))

TimeA_Text.on_submit(submitTimeA)

def submitTimeB(text):
    global TimeB
    if text == '':
        TimeB=np.nan
    else:
        TimeB=text   
    print('End Time={}'.format(TimeB))

TimeB_Text.on_submit(submitTimeB)

def PlotTimeAVG(event):
    global TimeA
    global TimeB
    global AVG
    
    if TimeA < TimeB:
      AVG=1
      update(event)
      AVG=0
    else:
      print('Invalid values for Time Range!')
        
TimeAVGButton.on_clicked(PlotTimeAVG)

def arrowclick(event):
    if event.key == 'right' and TimeSlider.val<len(Data[i]['time']):
        TimeSlider.set_val(TimeSlider.val+1)
    elif event.key == 'left' and TimeSlider.val>0:
        TimeSlider.set_val(TimeSlider.val-1)
    else:
        pass        

cid = fig.canvas.mpl_connect('key_press_event', arrowclick)

    
