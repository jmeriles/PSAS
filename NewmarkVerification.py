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
import matplotlib.pyplot as plt
import scipy.interpolate as sp

ModelSet=['Planar','line']
#Number and define NODE locations
NODE=[]
NODE.insert(0,[0,0])
NODE.insert(1,[0,1])


CON=[[1,2]]

#Turn selfweight on and off
selfweight=0

#Add nodal masses
Node_Mass=[12,12]

# DO NOT EDIT THIS
Model=Model(NODE,CON,Node_Mass,selfweight,ModelSet)

#Set the Boundary conditions at the nodes where there are boundaries.
#Others will have no BC

Model.BOUN[0]=[1,1,1];


#USER Define loading
#Model.LOAD[2]=[0,-100,0]

#DO NOT EDIT THIS
# Given Node Locations, connectivity, Element Type, BC, and Loading, initiates elements
el=Model.init_el()

# USER Define E A and I, releases, q0 and w here

#w is defined as positive outwards and negative into the beam
#Order of basic forces is [axial, torsion, Myi,Myj,Mzi,Mzj]

for i in range (Model.numel):
    el[i].Iz=1
    el[i].A=1
    el[i].E= 39.47841760435743/3
    el[i].rho=0

#el[0].REL=[np.inf,np.inf,np.inf,np.inf,np.inf,0];
#el[2].REL=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf];



#DO NOT EDIT
el=Model.elMat(el)
Model.Create_B(el)
Model.Create_K()
Model.Create_M()

M = Model.M
K = Model.K

masslessDOF = []
massDOF = []
for i in range(len(M)):
    if M[i,i] < 1E-5:
        masslessDOF.append(i)
    else:
        massDOF.append(i)
        
#Model.M = Model.M[0:2,0:2]
#Model.K = Model.K[0:2,0:2]
#Model.FDOF = Model.FDOF[0:2]

NewOrder = np.hstack([massDOF,masslessDOF])
K_cond = np.matrix([[K[i,j] for i in NewOrder] for j in NewOrder])
K_mass = K_cond[0:len(massDOF),0:len(massDOF)]
K_mass_nmass = K_cond[len(massDOF):,0:len(massDOF)]
K_nmass = K_cond[len(massDOF):,len(massDOF):]
K_cond = K_mass-K_mass_nmass.T @ np.linalg.inv(K_nmass) @ K_mass_nmass
M_cond = np.matrix([[M[i,j] for i in NewOrder] for j in NewOrder])
M_cond = M_cond[0:len(massDOF),0:len(massDOF)]
Model.M = M_cond
Model.K = K_cond
Model.FDOF = Model.FDOF[0:len(massDOF)]

Ground_Motion='El_Centro.txt'
DynOp=[.5,0,.001,[1,0,0]]
#DO NOT EDIT
el=Model.elMat(el)
An=Analyze()

#Static Analysis
#An.Displ(Model,el)

#Dynamic analysis
Model.Eig()
Model.Create_C("Rayleigh",[2,0,0])
GM=GMread(Ground_Motion)
#GM = np.array(GM)

dtold = .02
told = np.arange(0,dtold*len(GM),dtold)


#interpolate eq
GM_new=sp.interp1d(told,GM)

dt=.001
time = np.arange(0,dtold*(len(GM)-1),dt)
GM = GM_new(time)
#Eq_a = Eq_a/40
GM = np.hstack((GM))
time = np.arange(0,dt*(len(GM)),dt)


#P = Model.P
#P = np.tile(P,(1,1000))
[Model.Udy,Model.U_dotdy,Model.U_ddotdy,Model.Pbar]=An.Dyn_Newmark(Model,DynOp,[],GM)

# T = len(Model.Udy[0,:])*.01
# time = np.arange(0,T,.01)
plt.plot(time,Model.Udy[0,:])
# #plt.xlim([0,30])


# plot=PlotModel()
# plot.Struct(Model,el)
# #plot.DefShape(Model,el)
# animation = plot.AnimateDefShape(Model,el,1)



