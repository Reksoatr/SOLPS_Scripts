# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
from TOOLS import TANH
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button, CheckButtons

Shot = '12'
Attempt = ['101'] #,'102','103','104','105']
PS=['.','.','.','.','.','-']

NeuDen = SOLPSPLOT(Shot,Attempt,['Ne','NeuDen'],PlotScheme='x',EXP=False) #,AVG=True,PlotScheme=PS)
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
CoreBound = NeuDen.KW['CoreBound']
Rmax = 0.01
Rmin = -0.01
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']
RR,Rexp,Rstr = NeuDen.GetRadCoords('rrsep',[0,0])

f0 = JXA
p0 = [0,3.5e20,0.005,1e18,1e21]
x0 = []
xi = []
log = 2
axcolor = 'lightgoldenrodyellow'
FIT=0
#depth=0.02
#yfit = np.zeros(RR.shape[0],CoreBound[1]-CoreBound[0])

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(3,2,width_ratios=[5,3],height_ratios=[4,4,1],hspace=0.3)

axcontour = fig.add_subplot(gs[:,1])
NeuDen.Contour('NeuDen',LOG10=1,AX=axcontour, Markers=False)
axcontour.set_title('Neutral Density Contour')
l, = axcontour.plot(RadLoc.loc[:,f0,Attempt[-1]],VertLoc.loc[:,f0,Attempt[-1]],color='Red',linewidth=3)
axcontour.margins(x=0)

neudenprofile = fig.add_subplot(gs[0]) #plt.axes([0.25, 0.2, 0.4, 0.6], facecolor=axcolor)
NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=f0)  #,PlotScheme=['x'])
XLim = neudenprofile.get_xlim()
YLim = neudenprofile.get_ylim()

neprofile = fig.add_subplot(gs[2])
NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=f0)
neprofile.set_xlim(XLim)
NeYLim = neprofile.get_ylim()
neprofile.set_title('')

#ri = np.where(np.abs(RR.loc[:,f0,Attempt].values) > Rmin)[0][0]
#rf = np.where(np.abs(RR.loc[:,f0,Attempt].values) < Rmax)[0][-1]
#NDTrial = np.polyfit(RR.loc[ri:rf,f0,Attempt],np.log(NeuDen.PARAM['NeuDen'].loc[ri:rf,f0,Attempt]),1,full=True)

axslide = fig.add_subplot(gs[4], facecolor=axcolor) #plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sslide = Slider(axslide, 'Poloidal\nSurface', CoreBound[0], CoreBound[1]-1, valinit=f0, valfmt='%0.0f', valstep=1.0)
#axslide.set

def update(val):
    
    log=2
    
    if LogOn.get_status()[0] is True:
        log = 2
    elif LogOn.get_status()[0] is False:
        log = 0
        
    PolPos = sslide.val
    l.set_xdata(RadLoc.loc[:,PolPos,Attempt[-1]])
    l.set_ydata(VertLoc.loc[:,PolPos,Attempt[-1]])
    neudenprofile.clear()
    neprofile.clear()
    NeuDen.RadProf('NeuDen',LOG10=log,AX=neudenprofile,Markers=False,RADC='rrsep',JXA=PolPos) #,PlotScheme=['x'])
    NeuDen.RadProf('Ne',AX=neprofile,Markers=False,RADC='rrsep',JXA=PolPos)
    neudenprofile.set_xlim(XLim)
    neudenprofile.set_ylim(YLim)
    neprofile.set_xlim(XLim)
    neprofile.set_ylim(NeYLim)
    neprofile.set_title('')
    fig.canvas.draw_idle()


sslide.on_changed(update)

logonax = plt.axes([0.125, 0.025, 0.1, 0.06], facecolor=axcolor)
LogOn = CheckButtons(logonax, ['Log 10 Scale'],[1])
LogOn.on_clicked(update)

resetax = plt.axes([0.25, 0.025, 0.1, 0.06])
Reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sslide.reset()
    
Reset.on_clicked(reset)

tanhfitax = plt.axes([0.375, 0.025, 0.1, 0.06])
TanhFit = Button(tanhfitax, 'Create Tanh Fit', color=axcolor, hovercolor='0.975')

def tanhfit(event):
    PolPos=sslide.val
    if FIT == 0:
        RR_SOLPS = RR.loc[:,PolPos,Attempt[0]].values
        Ne_SOLPS = NeuDen.PARAM['Ne'].loc[:,PolPos,Attempt[0]].values
        yfit=curve_fit(TANH, RR_SOLPS, Ne_SOLPS,p0)
        yparam = yfit[0]
        print('Poloidal Slice {:0.0f}: r0={:.3f}m, h={:.3e}m^-3, d={:.3f}m, b={:.3e}m^-3, m={:.3e}m^-4'.format(PolPos,*yparam))
        neprofile.plot(RR_SOLPS,TANH(RR_SOLPS,*yparam))
        x0.append(yparam[0]+yparam[2])
        xi.append(yparam[0]-yparam[2])
        neprofile.axvline(x0[-1])
        neudenprofile.axvline(x0[-1])
        neprofile.axvline(xi[-1])
        neudenprofile.axvline(xi[-1])
        
        '''
        for n, i in enumerate(range(CoreBound[0],CoreBound[1])):
            yfit[:,n]=curve_fit(TANH, RR.loc[:,i,Attempt[0]].values, NeuDen.PARAM['Ne'].loc[:,i,Attempt[0]].values)
        '''
        
TanhFit.on_clicked(tanhfit)

expfitax = plt.axes([0.5, 0.025, 0.1, 0.06])
ExpFit = Button(expfitax, 'Create Exponential Fit', color=axcolor, hovercolor='0.975')

def expfit(event):
    PolPos=sslide.val
    if FIT == 0:
        xr = x0[-1]
        xri = xi[-1]
        RR_SOLPS = RR.loc[:,PolPos,Attempt[0]].values
        RR_i = np.where((RR_SOLPS>(xri)) & (RR_SOLPS<(xr)))[0]
        RR_SOLPS=RR_SOLPS[RR_i]
        NeuDen_SOLPS = NeuDen.PARAM['NeuDen'].loc[RR_i,PolPos,Attempt[0]].values
        efold=np.polyfit(RR_SOLPS,np.log(NeuDen_SOLPS),1,full=True)
        eparam = efold[0]
        print('Poloidal Slice {:0.0f}: e-folding length={:.3f}mm'.format(PolPos,1000/eparam[0]))
        print('Exponential fit from r-r_sep={:.3f}m to r-r_sep={:.3f}m'.format(RR_SOLPS[0],RR_SOLPS[-1]))
        NeuDenFit = np.exp(eparam[1]) * np.exp(eparam[0]*RR_SOLPS)
        neudenprofile.plot(RR_SOLPS, NeuDenFit)
        #neprofile.axvline(RR_SOLPS[0])
        #neudenprofile.axvline(RR_SOLPS[0])
        
ExpFit.on_clicked(expfit)

plt.show()