# -*- coding: utf-8 -*-
"""
B2.Transport.Inputfile Point Generator

Created on Wed Jul  4 11:48:11 2018

Version 1.2 - STABLE - Modified Tue Sep 18 12:25:00 2018

Edited 7/31/18 17:20:00
- Version 1 STABLE
Edited 8/23/18 15:06:00
- Extract function rewritten to use regex, outputs CoeffID to Combobox
- Version 1.1 STABLE
Edited 9/18/18 12:25:00
- Modified grid spacing of textbox and graph plot to be more compact
- Version 1.2 STABLE
Last Edited 11/8/2021 17:00:00
- Added in profile multiplication feature
- Version 1.3 STABLE

@author: Richard Reksoatmodjo
"""
import math
import re
import matplotlib

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2TkAgg
from matplotlib.backend_bases import MouseEvent
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import scrolledtext
from tkinter import *
from tkinter.ttk import *

class B2PtGen(tk.Tk): 
    
    def __init__(self):
        
        root = Tk()
        
        root.title("B2.Transport.Inputfile Point Generator")
        root.geometry("1000x800")
        
        app = Frame(root)
        app.grid()
        self._Xmin = -0.1
        self._Xmax = 0.1
        self._Ymin = 0
        self._Ymax = 1
        self._figure, self._axes, self._line = None, None, None
        self._newpoints = {}
        self._dragging_point = None
        self._points = {}
        
        Instr = Label(app,text = "Enter ranges of X and Y axes in corresponding text boxes below, then click 'Update Axes'. \nUse plot area to draw shape of Inputfile Transport Coefficient Profile \n(Left Mouse Button to create points, Right Mouse Button to delete points). \nThen click 'Generate Point List' to create Inputfile Point List", justify=CENTER)
        Instr.grid(column=0,row=0,columnspan=4)
        
        TextInstr = Label(app,text = "To Import data points, Copy and Paste \nONE set of transport coefficient points. \nInclude 'ndata' header line. Then click 'Extract Data Points'", justify=CENTER)
        TextInstr.grid(column=4,row=0)
        
        Lbl1 = Label(app, text = "X min")
        Lbl1.grid(column=0,row=1)
        
        XminT = Entry(app, width=10)
        XminT.insert(0,'-0.1')
        XminT.grid(column=0,row=2)
        
        Lbl2 = Label(app, text = "X max")
        Lbl2.grid(column=1,row=1)
        
        XmaxT = Entry(app, width=10)
        XmaxT.insert(0,'0.1')
        XmaxT.grid(column=1,row=2)
        
        Lbl3 = Label(app, text = "Y min")
        Lbl3.grid(column=0,row=3)
        
        YminT = Entry(app, width=10)
        YminT.insert(0,'0')
        YminT.grid(column=0,row=4)
        
        Lbl4 = Label(app, text = "Y max")
        Lbl4.grid(column=1,row=3)
        
        YmaxT = Entry(app, width=10)
        YmaxT.insert(0,'1')
        YmaxT.grid(column=1,row=4)
        
        def clicked():
            if float(XminT.get()) < float(XmaxT.get()) and float(YminT.get()) < float(YmaxT.get()): 
                self._Xmin = float(XminT.get())
                self._Xmax = float(XmaxT.get())
                self._Ymin = float(YminT.get())
                self._Ymax = float(YmaxT.get())
                Lbl1.configure(text = "Xmin = " + XminT.get())
                Lbl2.configure(text = "Xmax = " + XmaxT.get())
                Lbl3.configure(text = "Ymin = " + YminT.get())
                Lbl4.configure(text = "Ymax = " + YmaxT.get())
                
                self._init_plot(app)
                
            else:
                print("Input values do not form valid ranges")                  
        
        button1 = Button(app, command=clicked)
        button1.grid(column=2,row=1,columnspan=2)
        button1['text'] = "Update Axes"
        
        J = IntVar()
        J.set(1)
        
        AtomsRad = Radiobutton(app,text='Atoms',variable=J,value=0)
        IonsRad = Radiobutton(app,text='Ions',variable=J,value=1)
        AtomsRad.grid(column=2,row=4)
        IonsRad.grid(column=3,row=4)           
        
        Coeff = Combobox(app, width = 45, font=("Bold",12))
        Coeff['values'] = ('1 [dna0] particle density-driven diffusivity',
                           '2 [dpa0] particle pressure-driven diffusivity',
                           '3 [hcib] ion thermal anomalous diffusivity',
                           '4 [hce0] electron thermal anomalous diffusivity',
                           '5 [vla0x] Poloidal-component of the anomalous ”pinch” velocity',
                           '6 [vla0y] Radial-component of the anomalous ”pinch” velocity',
                           '7 [vsa0] anomalous viscosity',
                           '8 [sig0] anomalous radial electrical conductivity',
                           '9 [alf0] anomalous radial thermo-electric coefficient')
        Coeff.current(0)
        Coeff.grid(row=6, columnspan=4)

        def Generate():
            #print(self._points)
            n = len(self._points)
            m = 0
            i = str(Coeff.get())[0]
            j = J.get()
            r = sorted(self._points.items())
            print(' ndata(1, {0}, {1})= {2},'.format(i,j,n))
            for m in range(n):
                print(' tdata(1, {0}, {1}, {2})= {3}, tdata(2, {0}, {1}, {2})= {4},'.format(m+1,i,j,round(r[m][0],5),round(r[m][1],5)))
                
        button2 = Button(app, command=Generate)
        button2.grid(column=2,row=2,columnspan=2)
        button2['text'] = "Generate Point List"
        
        def Extract():
            TT = Textbox.get("1.0","end-1c")
            dataList = TT.split('\n')
            ndata = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", dataList[0])
            CoeffID = ndata[1]
            PtNo = ndata[3]
            XList = []
            YList = []
            for mm in range(int(PtNo)):
                XList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+1])[4]))
                YList.append(float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",dataList[mm+1])[9]))
            self._newpoints = dict(zip(XList,YList))
            Coeff.current(int(CoeffID)-1)
            print(self._newpoints)
            self._update_plot()
        
        button3 = Button(app, command=Extract)
        button3.grid(column=2,row=3,columnspan=2)
        button3['text'] = "Extract Data Points"
        
        Textbox = scrolledtext.ScrolledText(app, width=50,height=20)
        Textbox.grid(row=1,column=4,rowspan=6,sticky=N)
        
        def Multiply():
           M=float(MultiplierT.get())
           r = sorted(self._points.items())
           self._newpoints = dict(zip([x[0] for x in r],[M*y[1] for y in r]))
           self._update_plot()
           print('Values multiplied!') 
        
        button4 = Button(app, command=Multiply)
        button4.grid(column=4,row=7)
        button4['text'] = "Multiply Values:"
        
        MultiplierT = Entry(app, width=10)
        MultiplierT.insert(0,'1')
        MultiplierT.grid(column=4,row=8)
          
        root.mainloop() 
        
    def _init_plot(self, app):
        if not self._figure:
            self._figure = plt.figure(num=1,frameon=True)
            canvas = FigureCanvasTkAgg(self._figure, app)
            canvas.draw()
            canvas.get_tk_widget().grid(column=0,row=7,columnspan=4,rowspan=4)
            
        if not self._axes:
            self._axes = plt.axes()
                
        self._axes.set_xlim(self._Xmin, self._Xmax)
        self._axes.set_xlabel('Radial Distance from Separatrix (along Outer Midplane) [m]')
        self._axes.set_ylabel('Normalized Coefficient Magnitude')
        self._axes.set_ylim(self._Ymin, self._Ymax)
        self._axes.grid(b=True,which="both")
        self._figure.sca(self._axes)

        self._figure.canvas.mpl_connect('button_press_event', self._on_click)
        self._figure.canvas.mpl_connect('button_release_event', self._on_release)
        self._figure.canvas.mpl_connect('motion_notify_event', self._on_motion)
            
        self._figure.canvas.draw()
        
    def _update_plot(self):
        if self._newpoints:
            self._points = self._newpoints
        if not self._points:
            return
        x, y = zip(*sorted(self._points.items()))
        # Add new plot
        if not self._line:
            self._line, = self._axes.plot(x, y, "b", marker="o", markersize=5)
        # Update current plot
        else:
            self._line.set_data(x, y)
        self._figure.canvas.draw()

    def _add_point(self, x, y=None):
        if isinstance(x, MouseEvent):
            x, y = float(x.xdata), float(x.ydata)
        self._points[x] = y
        return x, y

    def _remove_point(self, x, _):
        if x in self._points:
            self._points.pop(x)

    def _find_neighbor_point(self, event):
        u""" Find point around mouse position
        :rtype: ((int, int)|None)
        :return: (x, y) if there are any point around mouse else None
        """
        distance_threshold = 0.05*(self._Ymax - self._Ymin)
        nearest_point = None
        min_distance = math.sqrt((self._Xmax - self._Xmin)**2 + (self._Ymax - self._Ymin)**2)
        for x, y in self._points.items():
            distance = math.hypot(event.xdata - x, event.ydata - y)
            if distance < min_distance:
                min_distance = distance
                nearest_point = (x, y)
        if min_distance < distance_threshold:
            return nearest_point
        return None

    def _on_click(self, event):
        u""" callback method for mouse click event
        :type event: MouseEvent
        """
        # left click
        if event.button == 1 and event.inaxes in [self._axes]:
            point = self._find_neighbor_point(event)
            if point:
                self._dragging_point = point
                self._remove_point(*point)
            else:
                self._add_point(event)   
            self._update_plot()
        # right click
        elif event.button == 3 and event.inaxes in [self._axes]:
            point = self._find_neighbor_point(event)
            if point:
                self._remove_point(*point)
                self._update_plot()

    def _on_release(self, event):
        u""" callback method for mouse release event
        :type event: MouseEvent
        """
        if event.button == 1 and event.inaxes in [self._axes] and self._dragging_point:
            self._add_point(event)
            self._dragging_point = None
            self._update_plot()

    def _on_motion(self, event):
        u""" callback method for mouse motion event
        :type event: MouseEvent
        """
        if not self._dragging_point:
            return
        self._remove_point(*self._dragging_point)
        self._dragging_point = self._add_point(event)
        self._update_plot()
        
if __name__ == "__main__":
    B2PtGen()            
        
