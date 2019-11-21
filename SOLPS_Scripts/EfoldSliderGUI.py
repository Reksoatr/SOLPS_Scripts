# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:15:42 2019

@author: rmreksoatmodjo
"""

from VesselPlotterNew import SOLPSPLOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider, Button, CheckButtons

Shot = 'gas12'
Attempt = '20'

NeuDen = SOLPSPLOT(Shot,Attempt,'NeuDen')
JXA = NeuDen.KW['JXA']
JXI = NeuDen.KW['JXI']
RadLoc = NeuDen.RadCoords['RadLoc']
VertLoc = NeuDen.RadCoords['VertLoc']

f0 = JXA
log = 2
axcolor = 'lightgoldenrodyellow'

fig, ax = plt.subplots()
ax.set_frame_on(False)
ax.set_axis_off()
gs=gridspec.GridSpec(2,2,width_ratios=[5,3],height_ratios=[8,1])

axcontour = fig.add_subplot(gs[:,1])
NeuDen.Contour('NeuDen',LOG10=1,AX=axcontour, Markers=False)
axcontour.set_title('Neutral Density Contour')
l, = axcontour.plot(RadLoc.loc[:,f0,Attempt],VertLoc.loc[:,f0,Attempt],color='Red',linewidth=3)
axcontour.margins(x=0)

axprofile = fig.add_subplot(gs[0]) #plt.axes([0.25, 0.2, 0.4, 0.6], facecolor=axcolor)
NeuDen.RadProf('NeuDen',LOG10=log,AX=axprofile,Markers=False,RADC='rrsep',JXA=f0,PlotScheme=['x'])

axslide = fig.add_subplot(gs[2], facecolor=axcolor) #plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sslide = Slider(axslide, 'Poloidal Surface', 1, 96, valinit=f0, valfmt='%0.0f', valstep=1.0)
#axslide.set

def update(val):
    '''
    if radio.get_status() is True:
        log = 2
    elif radio.get_status() is False:
        log = 0
    '''
    PolPos = sslide.val
    l.set_xdata(RadLoc.loc[:,PolPos,Attempt])
    l.set_ydata(VertLoc.loc[:,PolPos,Attempt])
    axprofile.clear()
    NeuDen.RadProf('NeuDen',LOG10=log,AX=axprofile,Markers=False,RADC='rrsep',JXA=PolPos,PlotScheme=['x'])
    fig.canvas.draw_idle()


sslide.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sslide.reset()
    
button.on_clicked(reset)

rax = plt.axes([0.025, 0.025, 0.15, 0.15], facecolor=axcolor)
radio = CheckButtons(rax, ['LOG10'],[1])
radio.on_clicked(update)


plt.show()