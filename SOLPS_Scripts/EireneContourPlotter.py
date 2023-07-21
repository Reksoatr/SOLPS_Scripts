# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:34:21 2020

@author: Richard

V1.0 - Completed May 18, 2020

V2.0 - Updated June 17, 2023
"""

from TOOLS import SET_WDIR
import re
import linecache
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.path import Path
import matplotlib.gridspec as gridspec
from matplotlib.widgets import TextBox, Button, Slider, RadioButtons
from scipy.interpolate import griddata
from PARAMDICT import EireneDict
#from SOLPS_Plotter import SOLPSPLOT

Shot='1100308004'
Device='cmod'
Attempt='14Rf0.7_split2' # 14Rf0.7 for 1100308004, 18Rf0.6 for 1080416025, 24Rf2.0 for 1100305023 

MeshID='010'  # 026 used for Shot025, 020 used for Shot012, 001 for d3d, 
              # 025/010 for 1100308004, 024/009 for 1100305023, 027/011 for 1080416025
LOG=True
Pressure=True
Param='PDENA'
F_157='0'
B2=True
B2_Param='Ne'

TeVAC=0.02 #Hard-wired Electron Temperature in EIRENE vacuum cells
NeVAC=1e8 #Hard-wired Electron Density in EIRENE vacuum cells
ED2P=(2/3)*1.20173129e-18 #Conversion from Energy Density to Pressure (eV*m^-3 to mTorr)

P0=np.empty((2))
P0.fill(np.nan)  #[1.65,-0.65])
P1=np.empty((2))
P1.fill(np.nan)  #[2,-1.1])

#button1=plt.imread('Icons/rounded-rectangle-button.png')

BASEDRT, TOPDRT = SET_WDIR('solps-iter/runs/','')
if 'd3d' in Shot:
    BASEDRT = '{}{}'.format(BASEDRT,Shot)
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))
else:
    BASEDRT = '{}cmod/{}home'.format(BASEDRT, Shot)
    DRT = '{}/Attempt{}/EirOutput'.format(BASEDRT, str(Attempt))

tz=np.loadtxt('{}/TriangVertLoc{}'.format(DRT,Attempt),usecols = (2))

tr=np.loadtxt('{}/TriangRadLoc{}'.format(DRT,Attempt),usecols = (2))

try:
    WallFile=np.loadtxt('{}/mesh.extra'.format(BASEDRT))
    WF=True
except:
    print('mesh.extra file not found! Using vvfile.ogr instead') 
    WF=False
    VVFILE = np.loadtxt('{}/vvfile.ogr'.format(BASEDRT))

if Param == 'fort.157':
    file='{}/{}'.format(DRT,Param)
    Parameter=EireneDict[F_157]
    
    with open(file,'r') as fp:
        for i,L in enumerate(fp):
            x=re.search(Parameter['Header'], L)
            if x is not None:
                print(L)
                line_start=i+5
                if F_157=='0':
                    line_start+=1
                
    row_num=int(linecache.getline(file,line_start).split()[0])-1
    
    Data=(100*100*100*2.052E-17/(4*np.pi))*np.genfromtxt(file,skip_header=line_start,max_rows=row_num,usecols=(2))               
    
else:
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
    if Param[-1]=='M':
        Parameter['Label']=r'Molecular D2 Pressure $(mTorr)$'
    elif Param[-1]=='A':
        Parameter['Label']=r'Atomic D Pressure $(mTorr)$'    
    
if LOG:
    Data=np.ma.log10(Data)
    Data=Data.filled(np.floor(Data.min()))
    LogTxt=r'$Log_{10}$ of '
else:
    LogTxt=''

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
Profile.set_xlabel('Distance along Chord from P0 (m)')
Profile.set_ylabel(Parameter['Label'])

P0_point, =Contour.plot(np.nan,np.nan,'rx')
P1_point, =Contour.plot(np.nan,np.nan,'gx')
Chord1, =Contour.plot(np.nan,np.nan,'b-')
ChordA, =Contour.plot(np.nan,np.nan,'k-')
ChordB, =Contour.plot(np.nan,np.nan,'k-')
ChordXY, =Contour.plot(np.nan,np.nan,'r.')
Prof1, =Profile.plot(np.nan,np.nan,'k.')
Prof2, =Profile.plot(np.nan,np.nan,'r-')

#resetax = plt.axes([0.125, 0.025, 0.05, 0.05])
#Reset = Button(resetax, 'Reset', hovercolor='0.975') 

textP0Xax = plt.axes([0.11, 0.3, 0.05, 0.05])
P0X_Text = TextBox(textP0Xax, r'$R_{P0}$ (m)', hovercolor='0.9')

textP0Yax = plt.axes([0.22, 0.3, 0.05, 0.05])
P0Y_Text = TextBox(textP0Yax, r'$Z_{P0}$ (m)', hovercolor='0.9')

buttonP0ax = plt.axes([0.28, 0.3, 0.05, 0.05])
P0Button = Button(buttonP0ax, 'Set P0')#,image=button1)

clearP0ax = plt.axes([0.34, 0.3, 0.07, 0.05])
P0Clear = Button(clearP0ax, 'Clear P0')#,image=button1)

textP1Xax = plt.axes([0.11, 0.2, 0.05, 0.05])
P1X_Text = TextBox(textP1Xax, r'$R_{P1}$ (m)', hovercolor='0.9')

textP1Yax = plt.axes([0.22, 0.2, 0.05, 0.05])
P1Y_Text = TextBox(textP1Yax, r'$Z_{P1}$ (m)', hovercolor='0.9')

buttonP1ax = plt.axes([0.28, 0.2, 0.05, 0.05])
P1Button = Button(buttonP1ax, 'Set P1')#,image=button1)

clearP1ax = plt.axes([0.34, 0.2, 0.07, 0.05])
P1Clear = Button(clearP1ax, 'Clear P1')#,image=button1)

radioax = plt.axes([0.42, 0.2, 0.07, 0.15])
Radio = RadioButtons(radioax, ('nearest', 'linear', 'cubic'),active=1)

plotax = plt.axes([0.5, 0.2, 0.05, 0.15])
PlotChord = Button(plotax, 'Plot\nChord')#,image=button1)

axslide = plt.axes([0.15, 0.1, 0.4, 0.03])
sslide = Slider(axslide, 'Threshhold', 0.001, 0.025, valinit=0.010, valfmt='%0.3f', valstep=0.001)

if WF:
    for i in range(len(WallFile[:,0])):
        Contour.plot((WallFile[i,0],WallFile[i,2]),(WallFile[i,1],WallFile[i,3]),'k-',linewidth=3.0)
else:
    Contour.plot(VVFILE[:,0]/1000,VVFILE[:,1]/1000,'k-',linewidth=3.0)

IM=Contour.tripcolor(TP,Data)
Contour.set_aspect('equal')
Contour.grid()
Contour.set_xlabel('Radial Position R (m)')
Contour.set_ylabel('Vertical Position Z (m)')
Contour.set_title('{}{} \n Shot{} Attempt{}'.format(LogTxt,Parameter['Label'],Shot,Attempt))
plt.colorbar(IM,ax=Contour)

def ClearP0(event): 
    P0.fill(np.nan)
    P0X_Text.set_val('')
    P0Y_Text.set_val('')
    P0_point.set_data(P0)
    
P0Clear.on_clicked(ClearP0)

def ClearP1(event): 
    P1.fill(np.nan)
    P1X_Text.set_val('')
    P1Y_Text.set_val('')
    P1_point.set_data(P1)
    
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
        #print(event.xdata,event.ydata)
        P0[0]=round(event.xdata,3)
        P0[1]=round(event.ydata,3)

def pointclickP1(event):
    if event.button==1 and event.inaxes == Contour:
        #print(event.xdata,event.ydata)
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
    P0_point.set_data([P0[0]],[P0[1]])
    EirFig.canvas.draw()

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
    P1_point.set_data([P1[0]],[P1[1]])
    EirFig.canvas.draw()

P1Button.on_clicked(setP1)

def chordplot(event):    
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
    Band=np.ma.array(Data,mask=~Mask).compressed()
    
    if LOG:
        Avg_Band = 10**np.mean(Band)
    else:
        Avg_Band = np.mean(Band)
    
    print('Average {}: {:.3}'.format(Parameter['Label'],Avg_Band ))
    
    if LOG and Pressure:
        Prof1.set_data(Chord,10**Band)
    else:
        Prof1.set_data(Chord,Band)

    Xline=np.linspace(P0[0],P1[0],num=100)
    Yline=np.linspace(P0[1],P1[1],num=100)
    
    IS=Radio.value_selected
    
    print('Interpolation setting currently set to {}'.format(IS))
    
    Vals=griddata((tr,tz),Data,(Xline,Yline),method=IS)
    Dist=np.sqrt((Xline-P0[0])**2 + (Yline-P0[1])**2)
    
    if LOG and Pressure:
        Prof2.set_data(Dist,10**Vals)
    else:
        Prof2.set_data(Dist,Vals)
    
    Profile.relim()
    Profile.autoscale_view(True,True,True)
    Profile.set_title(r'Average {}: {:.3}'.format(Parameter['Label'],Avg_Band))
    Profile.legend(['Cell Data','Interpolated\nProfile ({})'.format(IS)])
    
    EirFig.canvas.draw()
    
    Output={'P0':P0,'Chord':(ChordX,ChordY),'Data':Band,'InterpChord':Dist,'InterpData':Vals}
    
    return Output   # To access Output dict in workspace, just call A=chordplot(0)
    
PlotChord.on_clicked(chordplot)

def arrowclick(event):
    if event.key == 'right' and round(sslide.val,3)<0.025:
        sslide.set_val(sslide.val+0.001)
    elif event.key == 'left' and round(sslide.val,3)>0.001:
        sslide.set_val(sslide.val-0.001)
    else:
        pass        

cid2 = EirFig.canvas.mpl_connect('key_press_event', arrowclick)

plt.show()
 