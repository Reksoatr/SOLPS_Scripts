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
from SOLPS_Plotter import SOLPSPLOT
import numpy as np
import pickle as pkl
import sys

BASEDRT, TOPDRT=SET_WDIR('','')

Device='cmod'

Shots=['1100308004','1080416025','1100305023']#['1080416025','1101014029']#['1100308004','1101014006']#['1100305023','1101014019']#['1120917011']#['1160718024','1160718025']#['1160718012','1160718013','1160718023']#   

Rad='psin'#'rmid'#

Filter = True

SOLPS=False
Attempts=[]#[15N]
KW={'JXA':40}
'''
KW={'JXA' : 27,
    'JXI' : 38,
    'SEP' : 13,
    'XDIM' : 66,
    'YDIM' : 26,
    'CoreBound' : [16,47]}
'''
        
Time0=70
AVG=0
PsinOffset=[-0.015,0.03,0]#[0.03,-0.013]#[-0.015,-0.02]
RmidOffset=[0,0,0]#[-0.01,-0.015]
TimeA=np.nan
TimeB=np.nan

CoreTS=True
CoreTime=1.0
LCFS_Te=[90,83,93]#[83,98]#[90,90]#[93,110]#

XLIM=[0.60,1.05]#[0.45,1.06]
NeLIM=[]#[0,7.5e19]
TeLIM=[]#[0,1500]

Colors=['red','green','blue']

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(2,1,height_ratios=[1,1],hspace=0.0)

ne_profile = fig.add_subplot(gs[0,0])
te_profile = fig.add_subplot(gs[1,0],sharex=ne_profile)
#time_slide = fig.add_subplot(gs[2,0])

if SOLPS:
    Sim=SOLPSPLOT(Shots[0],Attempts,ROOTSHOT='',EXP=False,Markers=False,**KW)
    Sim.RadProf('Ne',AX=ne_profile)
    Sim.RadProf('Te',AX=te_profile)

N = len(Shots)
Data={}
NeLine={}
TeLine={}
Coords={'rmid':'R (m)','psin':'$\psi_n$'}

if CoreTS:
    
    for n, i in enumerate(Shots):
        with open('{}gfileProcessing/cmod_files/{}_CORE.pkl'.format(TOPDRT,i),'rb') as f:
            u = pkl._Unpickler(f)
            u.encoding = 'latin1'
            CTS = u.load()
            #CTS=pkl.load(f)
        
        Data[i]={}
        Data[i]['ne']=np.array([CTS[0],CTS[1]*1e20,CTS[2]*1e20])
        Data[i]['te']=np.array([CTS[3],CTS[4]*1e3,CTS[5]*1e3])
        
        Data[i]['ne']=Data[i]['ne'][:,Data[i]['ne'][0,:].argsort()]
        Data[i]['te']=Data[i]['te'][:,Data[i]['te'][0,:].argsort()]            
        
        if LCFS_Te!=0:
            index=(np.abs(Data[i]['te'][1]-LCFS_Te[n])).argmin()
            PsinOffset[n]=1-Data[i]['te'][0][index]
            print('Psin offset for shot {} is {}'.format(i,PsinOffset[n]))
        
        if Filter:
            for v in ['ne','te']:
                for t in range(1,len(Data[i][v][1])):
                    if (Data[i][v][1][t]-Data[i][v][2][t]) > 2*(Data[i][v][1][t-1]):
                        Data[i][v][0][t]=Data[i][v][1][t]=Data[i][v][2][t]=np.nan
                        
        Time0=0
        Data[i]['time']=[CoreTime]
        CoreNe = ne_profile.errorbar(Data[i]['ne'][0]+PsinOffset[n],Data[i]['ne'][1],yerr=Data[i]['ne'][2],linestyle=':',capsize=5,marker='.',color=Colors[n])
        CoreTe = te_profile.errorbar(Data[i]['te'][0]+PsinOffset[n],Data[i]['te'][1],yerr=Data[i]['te'][2],linestyle=':',capsize=5,marker='.',color=Colors[n])
        ne_profile.fill_between(Data[i]['ne'][0]+PsinOffset[n],Data[i]['ne'][1]+Data[i]['ne'][2],Data[i]['ne'][1]-Data[i]['ne'][2],color=Colors[n],alpha=0.25)
        te_profile.fill_between(Data[i]['te'][0]+PsinOffset[n],Data[i]['te'][1]+Data[i]['te'][2],Data[i]['te'][1]-Data[i]['te'][2],color=Colors[n],alpha=0.25)

    textTimeAax = plt.axes([0, 0, 0, 0])
    TimeA_Text = TextBox(textTimeAax,'')

    textTimeBax = plt.axes([0, 0, 0, 0])
    TimeB_Text = TextBox(textTimeBax,'')

    buttonTimeAVGax = plt.axes([0, 0, 0, 0])
    TimeAVGButton = Button(buttonTimeAVGax,'')

    slideTimeax = plt.axes([0, 0, 0, 0])
    TimeSlider = Slider(slideTimeax,'',0,1)
    TimeSlider.valtext.set_text('') 
    
else:
    
    for n,i in enumerate(Shots):
        Data[i]=loadmat('{}gfileProcessing/{}_files/{}.mat'.format(TOPDRT,Device,i))
        if Rad not in Data[i].keys():
            sys.exit('Error! Coordinate {} does not exist! Aborting!'.format(Rad))
        if PsinOffset!=0 and Rad=='psin':
            Data[i]['psin']=Data[i]['psin']+PsinOffset[n]
        if RmidOffset!=0 and Rad=='rmid':
            Data[i]['rmid']=Data[i]['rmid']+RmidOffset[n]
        
        Data[i]['time']=Data[i]['time'].flatten()
        NeLine[i] = ne_profile.errorbar(Data[i][Rad][:,Time0],Data[i]['ne'][:,Time0],yerr=Data[i]['nerr'][:,Time0],marker='.',linestyle=':',color=Colors[n])
        TeLine[i] = te_profile.errorbar(Data[i][Rad][:,Time0],Data[i]['te'][:,Time0],yerr=Data[i]['terr'][:,Time0],marker='.',linestyle=':',color=Colors[n]) 
    
    textTimeAax = plt.axes([0.15, 0.05, 0.03, 0.05])
    TimeA_Text = TextBox(textTimeAax, 'Initial\nTime', hovercolor='0.9')

    textTimeBax = plt.axes([0.2, 0.05, 0.03, 0.05])
    TimeB_Text = TextBox(textTimeBax, 'End\nTime', hovercolor='0.9')

    buttonTimeAVGax = plt.axes([0.25, 0.05, 0.11, 0.05])
    TimeAVGButton = Button(buttonTimeAVGax, 'Plot Time Averaged Profiles')

    slideTimeax = plt.axes([0.45, 0.05, 0.45, 0.05])
    TimeSlider = Slider(slideTimeax, 'Time', 0, len(Data[i]['time']), valinit=Time0, valfmt="%i", valstep=1)
    TimeSlider.valtext.set_text('{:0.3f}\nseconds'.format(Data[i]['time'][Time0]))
    
ne_profile.set_title('Thompson Scattering Profiles at {:0.3f} sec'.format(Data[i]['time'][Time0]))
ne_profile.set_ylabel(r'Electron Density $n_e\;(m^{-3})$')
ne_profile.legend([*Attempts,*Shots])
ne_profile.axhline(0.0,color='k')
ne_profile.grid()

te_profile.set_ylabel(r'Electron Temperature $T_e\;(eV)$')
te_profile.set_xlabel(Coords[Rad])
te_profile.legend([*Attempts,*Shots])
te_profile.axhline(0.0,color='k')
te_profile.grid()

def update(event):
   global TimeA
   global TimeB
   global AVG

   for i in Shots:
       if CoreTS:
           TitleText='{:0.3f}'.format(Data[i]['time'][0])
           continue
       elif AVG==1:
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
   ne_profile.set_title('Edge Thompson Scattering Profiles at {} seconds'.format(TitleText))    
   te_profile.relim()
   te_profile.autoscale_view()
   
   fig.canvas.draw_idle()
   
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

if XLIM != []:
    ne_profile.set_xlim(*XLIM)
if NeLIM != []:
    ne_profile.set_ylim(*NeLIM)
if TeLIM != []:
    te_profile.set_ylim(*TeLIM)    


cid = fig.canvas.mpl_connect('key_press_event', arrowclick)

    
