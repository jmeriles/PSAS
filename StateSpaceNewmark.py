#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:09:29 2022

@author: juanmeriles
"""
from IPython import get_ipython
import numpy as np
from Model import Model
from PlotModel import PlotModel
#get_ipython().magic('reset -sf')
from Analyze import Analyze
from GMread import GMread
import scipy.interpolate as sp
import matplotlib.pyplot as plt

ModelSet=['Planar','line']
#Number and define NODE locations
NODE=[]
NODE.insert(0,[0,0])
NODE.insert(1,[0,5])
#NODE.insert(2,[10,5])
#NODE.insert(3,[10,0])

CON = [[1,2]]
#CON=[[1,2],[2,3],[3,4]]

#Turn selfweight on and off
selfweight=0

#Add nodal masses
#Node_Mass=[10000,10000,10000,10000]
Node_Mass = [12,12]
damping = 0.1

# DO NOT EDIT THIS
Model=Model(NODE,CON,Node_Mass,selfweight,ModelSet)

#Set the Boundary conditions at the nodes where there are boundaries.
#Others will have no BC

Model.BOUN[0]=[1,1,1];
#Model.BOUN[3] = [1,1,0];


#USER Define loading
#Model.LOAD[2]=[5,0,0,0,0,0]

#DO NOT EDIT THIS
# Given Node Locations, connectivity, Element Type, BC, and Loading, initiates elements
el=Model.init_el()

# USER Define E A and I, releases, q0 and w here

#w is defined as positive outwards and negative into the beam
#Order of basic forces is [axial, torsion, Myi,Myj,Mzi,Mzj]

for i in range (Model.numel):
    el[i].Iz=10000
    el[i].A=1000000000000
    el[i].E=1
    el[i].rho=0

#el[1].Iz = 1000000000000000
#el[1].REL=[np.inf,np.inf,0,np.inf,np.inf,np.inf];

#el[0].w_y=10
el[0].hybrid = 1
#el[1].hybrid = 1
el[0].controlDOF = [0,0,0,0,1,0]


hystiff = 3.3
el[0].ke = np.array([[ 0,0,0,0,0,0],
                    [ 0,0,0,0,0,0],
                    [ 0,0,0,0,0,0],
                    [0,0,0,1.0e+6,0,0],
                    [ 0,0,0,0, hystiff, 0],
                    [ 0,0,0,0,0,0]])

# el[0].ke = np.array([[0,0,0,0,0,0],
#                       [0,0,0,0,0,0],
#                       [0,0,0,0,0,0],
#                       [0,0,0,0,0,0],
#                       [0,0,0,0,0,0],
#                       [0,0,0,0,0,0]])



Ground_Motion='El_Centro.txt'
DynOp=[.25,.5,.001,[1,0,0]]
#DO NOT EDIT
el=Model.elMat(el)
Model.Create_B(el)
Model.Create_K()
Model.Create_M()
Model.Create_Loading(el)
Model.CheckHybridEl(el)


An=Analyze()


Model.Eig()
Model.Create_C("Rayleigh",[damping,0,0])
C=np.array(Model.C)
K=np.array(Model.K)
M=np.array(Model.M)

GM=GMread(Ground_Motion)
dtold = .02
told = np.arange(0,dtold*len(GM),dtold)


#interpolate eq
GM_new=sp.interp1d(told,GM)

dt=.001
time = np.arange(0,dtold*(len(GM)-1),dt)
GM = GM_new(time)
GM = GM/40
GM = np.hstack((GM))
time = np.arange(0,dt*(len(GM)),dt)

flag = 0
P = []
if len(P) == 0:
    P=np.zeros((len(Model.FDOF),GM.shape[0]))
    g=386.4#in/s^2, this is the standard but I should probably let people change it
    GM=g*GM
    LoadingLen = GM.shape[0]
    flag = 1
if len(GM) == 0:
    GM=np.zeros((1,P.shape[1]))
    LoadingLen = P.shape[1]
    Res = np.zeros((1,P.shape[0]))        
influence=[]


if Model.ModelSet[0]=='Planar':
    for i in range(len(Model.FDOF)):
        if (Model.FDOF[i]+3) % 3 == 0 and DynOp[3][0]==1:
            influence.append(1)
        elif (Model.FDOF[i]+3-1) % 3==0 and DynOp[3][1]==1:
            influence.append(1)
        else:
            influence.append(0)
            
else:
    for i in range(len(Model.FDOF)):
        if (Model.FDOF[i]+6) % 6 == 0 and DynOp[3][0]==1:
            influence.append(1)
        elif (Model.FDOF[i]+6-1) % 6==0 and DynOp[3][1]==1:
            influence.append(1)
        elif (Model.FDOF[i]+6-2) % 6==0 and DynOp[3][2]==1:
            influence.append(1)
        else:
            influence.append(0)

influence = np.array([influence])
#Set of influence vectors for each fixed dof

#E = np.transpose(influence)
u = np.zeros(len(Model.FDOF))
udot = np.zeros(len(Model.FDOF))
uddot = np.zeros(len(Model.FDOF))
u_pred = np.zeros(len(Model.FDOF))

alpha = 0
Beta = ((1-alpha)**2)/4
gamma = (1-2*alpha)/2


P0 = P[:,0]
P_next=P[:,0]-np.transpose(np.matmul(M,GM[0]*np.transpose(influence)))
m_eff = Model.M + (1+alpha)*dt*gamma*Model.C + dt**2*Beta*(1+alpha)*Model.K
u_pred = u+dt*udot+(dt**2)*(1-2*Beta)*uddot/2


U_hist = []
U_Kalman = []
U_hist.append(u)
#U_hist.append(u_pred)

P_resisting = Model.Kh@u_pred

P_old = P[:,0]
P_m = np.array([hystiff*u_pred[0],0])
P_tot = P_m+P_resisting


P_eff = (1-alpha)*P_next+alpha*P0-(1-alpha)*P_tot-alpha*P_old- ((1-alpha)*C*dt*(1-gamma)+alpha*(dt**2)*Beta*K)@uddot.T-C@udot.T
uddot_old = uddot
uddot = np.linalg.solve(m_eff,P_eff[0])

x = np.zeros(3*len(Model.FDOF))
#x[4] = dt**2*uddot[0]
#x[5] = dt**2*uddot[1]
print(x)
z = np.zeros(3*len(Model.FDOF))

gamma = 1/2

M_eff = M + C*gamma*dt
M_effinv = np.linalg.inv(M_eff)

n = len(Model.FDOF)
A = np.block([[np.eye(n),np.eye(n),.5*np.eye(n)],
              [-gamma*(dt**2)*M_effinv@K,np.eye(n)-gamma*dt*M_effinv@C-gamma*(dt**2)*M_effinv@K,(1-gamma)*(np.eye(n)-gamma*dt*M_effinv@C)-gamma*(dt**2)/2*M_effinv],
              [-(dt**2)*M_effinv@K,-dt*M_effinv@C-(dt**2)*M_effinv@K,-M_effinv@C*dt*(1-gamma)-(dt**2)/2*M_effinv]])

E = -np.hstack(np.array([np.zeros(n),(gamma*(dt**2)*M_effinv@M@influence.T).T[0],((dt**2)*M_effinv@M@influence.T).T[0]]))

xhist = np.zeros((len(GM)-1,6))
for i in range(1,len(GM)-1):
    x = A@x+E*GM[i-1]
    xhist[i] = x

plt.plot(xhist)



