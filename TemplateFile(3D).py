#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 19:37:16 2020

@author: jim53
"""

from IPython import get_ipython
import numpy as np
from Model import Model
from PlotModel import PlotModel
#get_ipython().magic('reset -sf')
from Analyze import Analyze
from GMread import GMread

ModelSet=['3D','line']
#Number and define NODE locations
NODE=[]
NODE.insert(0,[0,0,0])
NODE.insert(1,[2.5,0,-4.33])
NODE.insert(2,[5,0,0])
NODE.insert(3,[2.5,5,-2.165])

CON=[[1,4],[2,4],[3,4]]

#Turn selfweight on and off
selfweight=0

#Add nodal masses
Node_Mass=[1,2,3,4]

# DO NOT EDIT THIS
Model=Model(NODE,CON,Node_Mass,selfweight,ModelSet)

#Set the Boundary conditions at the nodes where there are boundaries.
#Others will have no BC

Model.BOUN[0]=[1,1,1,0,0,0];
Model.BOUN[1]=[1,1,1,0,0,0];
Model.BOUN[2]=[1,1,1,0,0,0];

#USER Define loading
Model.LOAD[3]=[0,100,0,0,0,0]

#DO NOT EDIT THIS
# Given Node Locations, connectivity, Element Type, BC, and Loading, initiates elements
el=Model.init_el()

# USER Define E A and I, releases, q0 and w here

#w is defined as positive outwards and negative into the beam
#Order of basic forces is [axial, torsion, Myi,Myj,Mzi,Mzj]

for i in range (Model.numel):
    el[i].Iy=20000
    el[i].Iz=10000
    el[i].A=1000
    el[i].E=1
    el[i].G=1/(2*1+.3)
    el[i].J=1
    el[i].rho=1/1000

#el[2].REL=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf];

el[0].w_y=10
el[2].hybrid = 0



Ground_Motion='El_Centro.txt'
DynOp=[.25,.5,.01,[1,0,0]]
#DO NOT EDIT
el=Model.elMat(el)
Model.Create_B(el)
Model.Create_K()
Model.Create_M()

plot=PlotModel()
plot.Struct(Model,el)

#Dynamic analysis
Model.Eig()
Model.Create_C("Rayleigh",[1,1,1])
An=Analyze()
GM=GMread(Ground_Motion)
[Model.Udy,Model.U_dotdy,Model.U_ddotdy,Pbar]=An.Dyn_Newmark(Model,DynOp,[],GM)

#Static Analysis
An.Displ(Model,el)


plot=PlotModel()
plot.Struct(Model,el)
plot.DefShape(Model,el)
animation = plot.AnimateDefShape(Model,el,100)
