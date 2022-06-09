# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 14:49:57 2022

@author: 18313
"""

from manim import *
import equilibrium as eq
import numpy as np
import matplotlib.pyplot as plt
from omfit_classes import omfit_solps

fig,ax=plt.subplots()
ax.axis('equal')

GF=eq.equilibrium(gfile='g1100305023.01000')
GF.plot_flux(axes=ax,colorbar=False,Nlines=25)


b2fgmtry = omfit_solps.OMFITsolps('b2fgmtry')
crx=b2fgmtry['crx'].reshape(4,38,98)[0,1:,1:]
cry=b2fgmtry['cry'].reshape(4,38,98)[0,1:,1:]

# Using crx[0], which is bottom left cell coordinate
# Starting index at 1,1 to remove left and lower boundary cells;
# total meshgrid size now 36x96 (37x97 points)

ax.plot(crx[:,:],cry[:,:],'k') #Perpendicular Surfaces
ax.plot(crx[18:,:].T,cry[18:,:].T,'b') #SOL Magnetic Surfaces
ax.plot(crx[:18,:24].T,cry[:18,:24].T,'r') #Inner PFR Magnetic Surfaces
ax.plot(crx[:18,24:72].T,cry[:18,24:72].T,'g') #Core Magnetic Surfaces
ax.plot(crx[:18,72:].T,cry[:18,72:].T,'r') #Outer PFR Magnetic Surfaces

ax.plot(crx[18,:].T,cry[18,:].T,color='black',linewidth=5)  #Separatrix
ax.plot(crx[-1,:].T,cry[-1,:].T,color='purple',linewidth=5) #SOL
ax.plot(crx[0,24:72].T,cry[0,24:72].T,color='yellow',linewidth=5) #Core
ax.plot(crx[0,:24].T,cry[0,:24].T,color='magenta',linewidth=5)  #Inner PFR
ax.plot(crx[0,72:].T,cry[0,72:].T,color='magenta',linewidth=5)  #Outer PFR
ax.plot(crx[:,0],cry[:,0],color='orange',linewidth=5) #Inner Target
ax.plot(crx[:,-1],cry[:,-1],color='orange',linewidth=5) #Outer Target
ax.plot(crx[:19,24],cry[:19,24],color='brown',linewidth=5,linestyle='dashed')   #Periodic Boundary
ax.plot(crx[:19,72],cry[:19,72],color='brown',linewidth=5,linestyle='dashed')   #Periodic Boundary
'''
#Border lines, if full 38x98 meshgrid was desired (39x99 points)
ax.plot(crx[1,:,-1],cry[1,:,-1],'k')
ax.plot(crx[1,:19,-2:].T,cry[1,:19,-2:].T,'r')
ax.plot(crx[1,19:,-2:].T,cry[1,19:,-2:].T,'b')

ax.plot(crx[2,-1,:].T,cry[2,-1,:].T,'b')
ax.plot(crx[2,-2:,:],cry[2,-2:,:],'k')

ax.plot(crx[3,-1,-2:].T,cry[3,-1,-2:].T,'b')
ax.plot(crx[3,-2:,-1],cry[3,-2:,-1],'k')
'''

'%%manim -qm SOLPSGridFull'

config.frame_height=1.0
config.frame_width=1.5
config.disable_caching=True

class SOLPSGridFull(Scene):
    def construct(self):
        x=96
        y=36

        M=np.meshgrid(range(x+1),range(y+1))
        xx=(M[0].flatten()-x/2)*0.0145
        yy=(M[1].flatten()-y/2)*0.0145

        x_trans=-0.65
        y_trans=0.12

        b2fgmtry = omfit_solps.OMFITsolps('b2fgmtry')
        crx=b2fgmtry['crx'].reshape(4,38,98)[0,1:,1:] + x_trans
        cry=b2fgmtry['cry'].reshape(4,38,98)[0,1:,1:] + y_trans

        background=ImageMobject("Cmod_Bfield")
        background.height=1.0
        background.set_x(0.1).set_y(y_trans)

        dots=VGroup().add(*[Dot([i,j,0],radius=0.001,fill_opacity=0.0) 
        for i,j in zip(crx.flatten(),cry.flatten())])

        dots2=VGroup().add(*[Dot([i,j,0],radius=0.001) 
        for i,j in zip(xx,yy)])

        plines=VGroup().set_z_index(0)        
        slines=VGroup().set_z_index(1)

        SOLlines=VGroup().set_z_index(1)
        PFRlines=VGroup().set_z_index(1)
        CORlines=VGroup().set_z_index(1)

        SOLbound=VGroup().set_z_index(2)
        SEPbound=VGroup().set_z_index(2)
        PFRbound=VGroup().set_z_index(2)
        CORbound=VGroup().set_z_index(2)

        ITARGbound=VGroup().set_z_index(3)
        OTARGbound=VGroup().set_z_index(3)
        PERbound=VGroup().set_z_index(3)

        PFRgap=VGroup()
        CORgap=VGroup()

        plines.add(*[Line(dots[i].get_center(),dots[1+i+x].get_center(),
        stroke_width=0.1) for i in range(y*(x+1))])

        for i in range(y*(x+1)):
            if np.mod(i,x+1)==0:
                ITARGbound.add(plines[i])
            elif np.mod(i+1,x+1)==0:
                OTARGbound.add(plines[i])
            elif np.mod(i-i//(x+1),x//4)==0 and np.mod(i-i//(x+1),x//2)!=0 and i < (y/2)*(x+1):
                PERbound.add(plines[i])

        for j in range(x*y//2):
            if j < x//4+x*(j//(x+1))-1:
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1))
                PFRlines.add(slines[j])
                if j < x+1:
                    PFRbound.add(slines[j])
            elif j == x//4+x*(j//(x+1))-1:
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1,color=PURPLE,stroke_opacity=0.0))
                if j < x+1:
                    slines[j].set_stroke(color=PINK,opacity=0.0).set_stroke_width(0.3)    
            elif j > x//4+x*(j//(x+1))-1 and j < 3*x//4+x*(j//(x+1))-1:
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1))
                CORlines.add(slines[j])
                if j < x+1:
                    CORbound.add(slines[j])
            elif j == 3*x//4+x*(j//(x+1))-1:
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1,color=TEAL,stroke_opacity=0.0))
                if j < x+1:
                    slines[j].set_stroke(color=TEAL_E,opacity=0.0).set_stroke_width(0.3)
            elif j > 3*x//4+x*(j//(x+1))-1:
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1))
                PFRlines.add(slines[j])
                if j < x:
                    PFRbound.add(slines[j])

        for j in range(x*y//2,x*y//2+x):    
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1))
                SEPbound.add(slines[j])

        for j in range(x*y//2+x,x*(y+1)):
                slines.add(Line(dots[j+j//x].get_center(),dots[1+j+j//x].get_center(),
                stroke_width=0.1))
                SOLlines.add(slines[j])
                if j >= x*y:
                    SOLbound.add(slines[j])

        SOLtext=Text("SOL").scale(0.03).next_to(SOLlines,RIGHT*0.1)
        PFRtext=Text("PFR").scale(0.03).move_to([-0.12,-0.35,0])
        CORtext=Text("Core").scale(0.03).align_to(CORlines).shift(UP*0.02)
        ITtext=Text("Inner\nTarget", line_spacing=1.5).scale(0.03).next_to(ITARGbound,LEFT*0.1)
        OTtext=Text("Outer\nTarget", line_spacing=1.5).scale(0.03).next_to(OTARGbound,RIGHT*0.1)
        IT2text=Text("Inner Target").scale(0.03)
        OT2text=Text("Outer Target").scale(0.03)
        SEPtext=Text("Separatrix").scale(0.03).next_to(SEPbound,LEFT*0.2)
        PERBCtext=Text("Periodic Boundary", line_spacing=1).scale(0.03)
        SOLBCtext=Text("SOL Boundary", line_spacing=1.5).scale(0.03)
        PFRBCtext=Text("PFR Boundary", line_spacing=1.5).scale(0.03)
        CORBCtext=Text("Core Boundary", line_spacing=1.5).scale(0.03)

        EXPLAINtext1=Text(
            "36 contours are traced\nalong surfaces of\nmagnetic equilibrium", line_spacing=2).scale(0.03).move_to([-0.5,0.34,0])

        EXPLAINtext2=Text(
            "Separatrix forms\nboundary between\nconfined Core region\nand unconfined\nscrape-off-layer (SOL)\nand private flux (PFR)\nregions", line_spacing=2).scale(0.03).move_to([0.5,0.26,0])

        EXPLAINtext3=Text(
            "Contours are split into\n96 poloidal sections,\nforming a curvi-\nlinear meshgrid of\n96 x 36 cells", line_spacing=2).scale(0.03).move_to([-0.5,0.3,0])

        EXPLAINtext4=Text(
            "Curvilinear meshgrid\nis transformed into\nrectangular domain\nby making a vertical\ncut from the Core to\nthe PFR through the\nX-point", line_spacing=2).scale(0.03).move_to([0.5,0.26,0])    

        self.play(FadeIn(background),run_time=3)
        self.wait
        self.play(background.animate.scale_to_fit_height(1.4))        
        self.play(FadeIn(EXPLAINtext1),Create(slines),run_time=3)
        self.wait
        self.play(FadeOut(background),run_time=3)
        self.wait
        self.play(FadeIn(EXPLAINtext2))
        self.wait(2)
        self.play(SEPbound.animate.set_stroke(color=PURE_RED,width=0.4),FadeIn(SEPtext))
        self.wait
        self.play(CORlines.animate.set_stroke(color=TEAL),FadeIn(CORtext))
        self.wait
        self.play(SOLlines.animate.set_stroke(color=BLUE),FadeIn(SOLtext))
        self.wait
        self.play(PFRlines.animate.set_stroke(color=PURPLE),FadeIn(PFRtext))
        self.wait
        self.play(ITARGbound.animate.set_stroke(color=PURE_GREEN,width=0.4),FadeIn(ITtext))
        self.wait
        self.play(OTARGbound.animate.set_stroke(color=PURE_GREEN,width=0.4),FadeIn(OTtext))
        self.wait(2)
        self.play(Transform(EXPLAINtext1,EXPLAINtext3))
        self.wait(2)
        self.play(Create(plines),run_time=3)
        self.wait(2)
        self.play(Transform(EXPLAINtext2,EXPLAINtext4))
        self.wait(2)
        self.play(PERbound.animate.set_stroke(color=YELLOW,width=0.6))        
        self.play(FadeOut(SEPtext,CORtext,SOLtext,PFRtext,ITtext,OTtext,
            EXPLAINtext1,EXPLAINtext2,EXPLAINtext3,
            EXPLAINtext4),run_time=2)
        self.play(CORbound.animate.set_stroke(color=TEAL_E,width=0.4),
            PFRbound.animate.set_stroke(color=PINK,width=0.4),
            SOLbound.animate.set_stroke(color=PURE_BLUE,width=0.4))

        for i,p in enumerate(plines):
            p.add_updater(lambda z,i=i,p=p: z.become(Line(dots[i].get_center(),
            dots[1+i+x].get_center(),color=p.get_color(),stroke_width=p.get_stroke_width())))

        for j,s in enumerate(slines):
            s.add_updater(lambda z,j=j,s=s: z.become(Line(dots[j+j//x].get_center(),
            dots[1+j+j//x].get_center(),color=s.get_color(),stroke_width=s.get_stroke_width(),
            stroke_opacity=s.get_stroke_opacity())))
            if j == x//4+x*(j//(x+1))-1:
                PFRgap.add(s)
            elif j == 3*x//4+x*(j//(x+1))-1:
                CORgap.add(s)            

        self.play(Transform(dots,dots2),run_time=5)
        self.wait(2)
        
        SOLtext.move_to([0,0.15,0])
        PFRtext.move_to([-0.53,-0.15,0])
        PFR2text=PFRtext.copy().move_to([0.53,-0.15,0])
        CORtext.move_to([0,-0.15,0])
        SEPtext.move_to([0,0.025,0])
        IT2text.move_to([-0.72,0,0]).rotate(PI/2)
        OT2text.move_to([0.72,0,0]).rotate(-PI/2)
        SOLBCtext.move_to([0,0.29,0])
        CORBCtext.move_to([0,-0.29,0])
        PFRBCtext.move_to([-0.53,-0.29,0])
        PFRBC2text=PFRBCtext.copy().move_to([0.53,-0.29,0])
        PERBCtext.move_to([0,-0.35,0])
        Arrow1=Arrow(start=PERBCtext.get_left(),end=CORbound[0].get_start())
        Arrow2=Arrow(start=PERBCtext.get_right(),end=CORbound[-1].get_end())
        
        self.play(FadeIn(SOLtext,PFRtext,PFR2text,CORtext,SEPtext,
            IT2text,OT2text,SOLBCtext,CORBCtext,
            PFRBCtext,PFRBC2text,PERBCtext,
            Arrow1,Arrow2),run_time=3)
        self.wait(2)
        
        
        