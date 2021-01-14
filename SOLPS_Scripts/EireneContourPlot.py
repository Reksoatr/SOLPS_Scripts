# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:34:21 2020

@author: Richard

V1.0 - Completed May 18, 2020
"""

from TOOLS import SET_WDIR
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.path import Path
import matplotlib.gridspec as gridspec
from matplotlib.widgets import TextBox, Button, Slider
from PARAMDICT import EireneDict

Shot='25'
Device='cmod'
Attempt='19N'
MeshID='026'
LOG=False
Pressure=True
Param='EDENA'

TeVAC=0.02 #Hard-wired Electron Temperature in EIRENE vacuum cells
NeVAC=1e8 #Hard-wired Electron Density in EIRENE vacuum cells
ED2P=(2/3)*1.20173129e-18 #Conversion from Energy Density to Pressure (eV*m^-3 to mTorr)

P0=np.empty((2))
P0.fill(np.nan)  #[1.65,-0.65])
P1=np.empty((2))  
P1.fill(np.nan)  #[2,-1.1])
#Thresh=0.05

button1=plt.imread('Icons/rounded-rectangle-button.png')

BASEDRT, TOPDRT = SET_WDIR('solps-iter/runs/','')
if 'd3d' in Shot:
    BASEDRT = '{}{}'.format(BASEDRT,Shot)
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))
else:
    BASEDRT = '{}cmod/0{}home'.format(BASEDRT, Shot[-2:])
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))

tz=np.loadtxt('{}/TriangVertLoc{}'.format(DRT,Attempt),usecols = (2))

tr=np.loadtxt('{}/TriangRadLoc{}'.format(DRT,Attempt),usecols = (2))

VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))

Parameter=EireneDict[Param]
Data=np.loadtxt('{}/{}{}'.format(DRT,Parameter['FileName'],Attempt),usecols = (2))

#NeuDen=np.loadtxt('{}/EirAtom{}'.format(DRT,Attempt),usecols = (2))
#MolDen=np.loadtxt('{}/EirMol{}'.format(DRT,Attempt),usecols = (2))
#TestIonDen=np.loadtxt('{}/EirTestIon{}'.format(DRT,Attempt),usecols = (2))
#NeuEng=np.loadtxt('{}/AtomEnergy{}'.format(DRT,Attempt),usecols = (2))
#MolEng=np.loadtxt('{}/MolEnergy{}'.format(DRT,Attempt),usecols = (2))
#TestIonEng=np.loadtxt('{}/TestIonEnergy{}'.format(DRT,Attempt),usecols = (2))

if 'EDEN' in Param and Pressure:
    Data=ED2P*Data
    Parameter['Label']=r'Pressure $(mTorr)$'
    
if LOG:
    Data=np.ma.log10(Data)
    Data=Data.filled(np.floor(Data.min()))

Nodes=np.fromfile('{}/{}.tria.{}.nodes'.format(BASEDRT,Device,MeshID),sep=' ') #Alternatively use fort.33
NN=int(Nodes[0])
XNodes=Nodes[1:NN+1]/100
YNodes=Nodes[NN+1:]/100

Triangles=np.loadtxt('{}/{}.tria.{}.cells'.format(BASEDRT,Device,MeshID),skiprows=1, usecols=(1,2,3)) #Alternatively use fort.34

TP=tri.Triangulation(XNodes,YNodes,triangles=(Triangles-1))

#####

EirFig, EirAx = plt.subplots()
EirFig.tight_layout()
EirAx.set_frame_on(False)
EirAx.set_axis_off()
gs=gridspec.GridSpec(2,2,width_ratios=[6,4],height_ratios=[4,3],hspace=0.05)
'''
Triangles=np.loadtxt('{}/triangles.dat'.format(BASEDRT))
Triangles=Triangles/100
NT=int(Triangles.shape[0]/4)
for n in range(0,NT):
    EirAx.plot(Triangles[n*4:(n*4)+3,0],Triangles[n*4:(n*4)+3,1],'k-')
'''
Contour = EirFig.add_subplot(gs[0:2,1])
Profile = EirFig.add_subplot(gs[0,0])
Profile.set_xlabel('Distance along Chord from Core Boundary (m)')
Profile.set_ylabel(Parameter['Label'])

P0_point, =Contour.plot(np.nan,np.nan,'rx')
P1_point, =Contour.plot(np.nan,np.nan,'gx')
Chord1, =Contour.plot(np.nan,np.nan,'b-')
ChordA, =Contour.plot(np.nan,np.nan,'k-')
ChordB, =Contour.plot(np.nan,np.nan,'k-')
ChordXY, =Contour.plot(np.nan,np.nan,'r.')
Prof1, =Profile.plot(np.nan,np.nan,'k.')

#resetax = plt.axes([0.125, 0.025, 0.05, 0.05])
#Reset = Button(resetax, 'Reset', hovercolor='0.975') 

textP0Xax = plt.axes([0.2, 0.3, 0.03, 0.05])
P0X_Text = TextBox(textP0Xax, r'$R_{P0}$ (m)', hovercolor='0.9')

textP0Yax = plt.axes([0.275, 0.3, 0.03, 0.05])
P0Y_Text = TextBox(textP0Yax, r'$Z_{P0}$ (m)', hovercolor='0.9')

buttonP0ax = plt.axes([0.32, 0.3, 0.05, 0.05])
P0Button = Button(buttonP0ax, 'Set P0',image=button1)

clearP0ax = plt.axes([0.38, 0.3, 0.07, 0.05])
P0Clear = Button(clearP0ax, 'Clear P0',image=button1)

textP1Xax = plt.axes([0.2, 0.2, 0.03, 0.05])
P1X_Text = TextBox(textP1Xax, r'$R_{P1}$ (m)', hovercolor='0.9')

textP1Yax = plt.axes([0.275, 0.2, 0.03, 0.05])
P1Y_Text = TextBox(textP1Yax, r'$Z_{P1}$ (m)', hovercolor='0.9')

buttonP1ax = plt.axes([0.32, 0.2, 0.05, 0.05])
P1Button = Button(buttonP1ax, 'Set P1',image=button1)

clearP1ax = plt.axes([0.38, 0.2, 0.07, 0.05])
P1Clear = Button(clearP1ax, 'Clear P1',image=button1)

plotax = plt.axes([0.45, 0.25, 0.1, 0.05])
PlotChord = Button(plotax, 'Plot Chord',image=button1)

axslide = plt.axes([0.15, 0.1, 0.4, 0.03])
sslide = Slider(axslide, 'Threshhold', 0.005, 0.1, valinit=0.05, valfmt='%0.3f', valstep=0.005)

Contour.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-')
IM=Contour.tripcolor(TP,Data)
Contour.set_aspect('equal')
Contour.grid()
Contour.set_xlabel('Radial Position R (m)')
Contour.set_ylabel('Vertical Position Z (m)')
Contour.set_title('{}, Shot{} Attempt{}'.format(Parameter['Label'],Shot,Attempt))
plt.colorbar(IM,ax=Contour)

def ClearP0(event): 
    P0.fill(np.nan)
    P0X_Text.set_val('')
    P0Y_Text.set_val('')
    P0_point.set_data(P0)

def ClearP1(event): 
    P1.fill(np.nan)
    P1X_Text.set_val('')
    P1Y_Text.set_val('')
    P1_point.set_data(P1)
    
P0Clear.on_clicked(ClearP0)
P1Clear.on_clicked(ClearP1)

def submitP0X(text):
    if text == '':
        P0[0]=np.nan
    else:
        P0[0]=text    
    print('P0X={}'.format(P0[0]))

P0X_Text.on_submit(submitP0X)

def submitP0Y(text):
    if text == '':
        P0[1]=np.nan
    else:
        P0[1]=text   
    print('P0Y={}'.format(P0[1]))

P0Y_Text.on_submit(submitP0Y)

def submitP1X(text):
    if text == '':
        P1[0]=np.nan
    else:
        P1[0]=text   
    print('P1X={}'.format(P1[0]))

P1X_Text.on_submit(submitP1X)
    
def submitP1Y(text):
    if text == '':
        P1[1]=np.nan
    else:
        P1[1]=text   
    print('P1Y={}'.format(P1[1]))

P1Y_Text.on_submit(submitP1Y)

def pointclickP0(event):
    if event.button==1 and event.inaxes == Contour:
        print(event.xdata,event.ydata)
        P0[0]=round(event.xdata,3)
        P0[1]=round(event.ydata,3)

def pointclickP1(event):
    if event.button==1 and event.inaxes == Contour:
        print(event.xdata,event.ydata)
        P1[0]=round(event.xdata,3)
        P1[1]=round(event.ydata,3)

def setP0(event):
    global P0
    global P0_point
    if np.all(np.isnan(P0)):
        print('Select Point!')
        cid=EirFig.canvas.mpl_connect('button_press_event',pointclickP0)
        EirFig.waitforbuttonpress(-1)
        EirFig.canvas.mpl_disconnect(cid)

    P0X_Text.set_val(P0[0])
    P0Y_Text.set_val(P0[1])
    P0_point.set_data(P0[0],P0[1])
    EirFig.canvas.draw()
    #print(P0[0],P0[1])

P0Button.on_clicked(setP0)

def setP1(event):
    global P1
    global P1_point
    if np.all(np.isnan(P1)):
        print('Select Point!')
        cid=EirFig.canvas.mpl_connect('button_press_event',pointclickP1)
        EirFig.waitforbuttonpress(-1)
        EirFig.canvas.mpl_disconnect(cid)

    P1X_Text.set_val(P1[0])
    P1Y_Text.set_val(P1[1])
    P1_point.set_data(P1[0],P1[1])
    EirFig.canvas.draw()
    #print(P1[0],P1[1])

P1Button.on_clicked(setP1)

def chordplot(event):
    global Prof1
    Thresh=round(sslide.val,3)
    print('Threshhold={}m'.format(Thresh))
    PP=P1-P0
    Theta=np.arctan(PP[1]/PP[0])
    displace=Thresh*np.array([-np.sin(Theta),np.cos(Theta)])
    P0A=P0+displace
    P0B=P0-displace
    P1A=P1+displace
    P1B=P1-displace
    
    Chord1.set_data([P0[0],P1[0]],[P0[1],P1[1]])
    ChordA.set_data([P0A[0],P1A[0]],[P0A[1],P1A[1]])
    ChordB.set_data([P0B[0],P1B[0]],[P0B[1],P1B[1]])
    
    Bounds=Path([P0A,P1A,P1B,P0B])
    Mask=Bounds.contains_points(np.vstack((tr,tz)).T)
    ChordX=np.ma.array(tr,mask=~Mask).compressed()
    ChordY=np.ma.array(tz,mask=~Mask).compressed()
    
    ChordXY.set_data(ChordX,ChordY)
    
    Chord=np.sqrt((ChordX-P0[0])**2 + (ChordY-P0[1])**2)
    Chord=Chord-Chord.min()
    Band=np.ma.array(Data,mask=~Mask).compressed()

    print('Distance along Chord:{}'.format(Chord))
    print('{}:{}'.format(Parameter['Label'],Band))
    print('Average {}: {}'.format(Parameter['Label'],np.mean(Band)))
    
    Prof1.set_data(Chord,Band)
    
    Profile.relim()
    Profile.autoscale_view(True,True,True)
    Profile.set_title(r'Average {}: {}'.format(Parameter['Label'],np.mean(Band)))
    
    EirFig.canvas.draw()
    
PlotChord.on_clicked(chordplot)

def arrowclick(event):
    if event.key == 'right' and round(sslide.val,3)<0.1:
        sslide.set_val(sslide.val+0.005)
    elif event.key == 'left' and round(sslide.val,3)>0.005:
        sslide.set_val(sslide.val-0.005)
    else:
        pass        

cid2 = EirFig.canvas.mpl_connect('key_press_event', arrowclick)

plt.show()

    
    