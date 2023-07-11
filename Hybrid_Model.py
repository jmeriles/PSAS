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
import scipy.interpolate as sp
import matplotlib.pyplot as plt
import scipy

ModelSet=['Planar','line']
#Number and define NODE locations
NODE=[]
NODE.insert(0,[0,0])
NODE.insert(1,[0,5])
NODE.insert(2,[-3,0])
NODE.insert(3,[3,0])

#CON = [[1,2]]
CON=[[1,2],[3,2],[4,2]]

#Turn selfweight on and off
selfweight=0

#Add nodal masses
#Node_Mass=[10000,10000,10000,10000]









Node_Mass = [5,5,5,5]
damping = 0.05









# DO NOT EDIT THIS
Model=Model(NODE,CON,Node_Mass,selfweight,ModelSet)

#Set the Boundary conditions at the nodes where there are boundaries.
#Others will have no BC

Model.BOUN[0]=[1,1,1];
Model.BOUN[2] = [1,1,0];
Model.BOUN[3] = [1,1,0];


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
                    [0,0,0,1.0e+8,0,0],
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


plot=PlotModel()
plot.Struct(Model,el)


An=Analyze()

#Static Analysis
# An.Displ(Model,el)
# plot = PlotModel()
# plot.Struct(Model,el)
# plot.DefShape(Model,el)




K=np.array(Model.K)
M=np.array(Model.M)


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

dt=.01
time = np.arange(0,dtold*(len(GM)-1),dt)
GM = GM_new(time)







GM = GM/5







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


u = np.zeros(len(Model.FDOF))
udot = np.zeros(len(Model.FDOF))
uddot = np.zeros(len(Model.FDOF))
u_pred = np.zeros(len(Model.FDOF))

u_true = np.zeros(len(Model.FDOF))
udot_true = np.zeros(len(Model.FDOF))
uddot_true = np.zeros(len(Model.FDOF))
u_pred_true = np.zeros(len(Model.FDOF))

alpha = 0
Beta = ((1-alpha)**2)/4
gamma = (1-2*alpha)/2


P0 = P[:,0]
P_next=P[:,0]-np.transpose(np.matmul(M,GM[0]*np.transpose(influence)))
#P_next=P[:,1]-np.transpose(np.matmul(M,GM[1]*np.transpose(influence)))
m_eff = Model.M + (1+alpha)*dt*gamma*Model.C + dt**2*Beta*(1+alpha)*Model.K



U_hist = []
U_hist_true = []
U_Kalman = []
U_hist.append(u)
U_hist_true.append(u)
#U_hist.append(u_pred)

P_resisting = Model.Kh@u_pred

P_old = P[:,0]
P_m = np.array([hystiff*u_pred[0],0])
P_tot = P_m+P_resisting
P_tot_true = P_m+P_resisting


P_eff = (1-alpha)*P_next+alpha*P0-(1-alpha)*P_tot-alpha*P_old- ((1-alpha)*C*dt*(1-gamma)+alpha*(dt**2)*Beta*K)@uddot.T-C@udot.T
uddot_old = uddot
uddot_old_true = uddot_true
uddot = np.linalg.solve(m_eff,P_eff[0])
uddot_true = np.linalg.solve(m_eff,P_eff[0])

#Initializing uddot here, everything else initialized as 0, not 100% sure this is right
#u = u_pred+(dt**2)*Beta*uddot
#udot = udot + dt*((1-gamma)*uddot_old+gamma*uddot)
#u_pred = u+dt*udot+(dt**2)*(1-2*Beta)*uddot/2

#U_hist.append(u)
T=dt*GM.shape[0]
t=np.linspace(0,T,GM.shape[0]-1)


P_hist = []
P_hist_true = []

#Kalman Filtering Steady State
# x = np.zeros(3*len(Model.FDOF))
# z = np.zeros(1)

# gammaK = 1/2

# M_eff = M + C*gammaK*dt
# M_effinv = np.linalg.inv(M_eff)

# n = len(Model.FDOF)
# A = np.block([[np.eye(n),np.eye(n),.5*np.eye(n)],
#               [-gammaK*(dt**2)*M_effinv@K,np.eye(n)-gammaK*dt*M_effinv@C-gammaK*(dt**2)*M_effinv@K,(1-gammaK)*(np.eye(n)-gammaK*dt*M_effinv@C)-gammaK*(dt**2)/2*M_effinv],
#               [-dt**2*M_effinv@K,-dt*M_effinv@C-(dt**2)*M_effinv@K,-M_effinv@C*dt*(1-gammaK)-(dt**2)/2*M_effinv]])

# E = -np.hstack(np.array([np.zeros(n),(gammaK*(dt**2)*M_effinv@M@influence.T).T[0],((dt**2)*M_effinv@M@influence.T).T[0]]))

# measureblock = np.ones(n)
# measureblock[0] = 1
# H = np.block([[K[0],np.zeros(n),np.zeros(n)]])


# N = len(x)
# V = np.zeros((N,N))
#V[0,0] = 0#.00000001
#V = 0.0000000000000001*np.eye(len(x))
# W = 33.33*np.eye(len(z))

# Pp=np.zeros((len(GM),N,N))
# Pm=np.zeros((len(GM),N,N))

# G=np.zeros((N,1,1))
# uKal = np.zeros(6)
# uKal[4] = dt**2*uddot[0]
# uKal[5] = dt**2*uddot[1]
# U_Kalman.append(uKal)

# Pp[0] = np.zeros(N)
# Pm[0] = np.zeros(N)

# Pinf = scipy.linalg.solve_discrete_are(A.T,H.T,V,W)
# Ginf = Pinf @ H.T @ np.linalg.inv(H @ Pinf @ H.T + W)

w1 = np.zeros(len(u))
w2 = 0

f = np.zeros((len(GM),2))
fprev = np.array([hystiff*u_pred[0],0])
f_true = np.zeros((len(GM),2))
fprev_true = np.array([hystiff*u_pred[0],0])







fy = 2








K2 = 0
du = 0

for i in range(len(GM)-1):

    
    du = u_pred[0]-u[0]
    if i ==0:
        f[i] = fprev+np.array([hystiff*du,0])
    else:
        f[i] = f[i-1] + np.array([hystiff*du,0])

    if f[i][0]>fy:
        f[i] = fprev+ np.array([K2*du,0])

    elif f[i][0]<-fy:
        f[i] = fprev+np.array([K2*du,0])
        
        
    du_true = u_pred_true[0]-u_true[0]
    if i ==0:
        f_true[i] = fprev_true+np.array([hystiff*du_true,0])
    else:
        f_true[i] = f_true[i-1] + np.array([hystiff*du_true,0])

    if f_true[i][0]>fy:
        f_true[i] = fprev_true+ np.array([K2*du_true,0])

    elif f_true[i][0]<-fy:
        f_true[i] = fprev_true+np.array([K2*du_true,0])
        
    fprev_true= f_true[i]
    fprev = f[i]
    
    Pi=np.copy(P_next)

    P_next=P[:,i+1]-np.transpose(np.matmul(M,GM[i+1]*np.transpose(influence)))
    w2 = (np.random.rand()-.5)*.05
    #P_m = np.array([hystiff*u_pred[0],0])+w2
    #P_m_true = np.array([hystiff*u_pred_true[0],0])
    P_m = np.array([f[i][0]+w2,0])
    #P_hist.append(P_m[0])
    P_m_true = np.array([f_true[i][0],0])

    P_old = np.copy(P_tot)
    P_old_true = np.copy(P_tot_true)
    P_resisting = Model.Kh@u_pred
    P_resisting_true = Model.Kh@u_pred_true


    P_tot = P_m+P_resisting
    P_tot_true = P_m_true+P_resisting_true
    
    P_hist.append(P_tot)
    P_hist_true.append(P_tot_true)
    
    
    P_eff = (1-alpha)*P_next+alpha*Pi-(1-alpha)*P_tot-alpha*P_old- ((1-alpha)*C*dt*(1-gamma)+alpha*(dt**2)*Beta*K)@uddot.T-C@udot.T
    P_eff_true = (1-alpha)*P_next+alpha*Pi-(1-alpha)*P_tot_true-alpha*P_old_true- ((1-alpha)*C*dt*(1-gamma)+alpha*(dt**2)*Beta*K)@uddot_true.T-C@udot_true.T

    z = P_tot
    #uKal = (np.eye(N)-Ginf@H)@A@uKal+(np.eye(N)-Ginf@H)@E*GM[i]+(Ginf.T*P_tot[0])[0]
    #U_Kalman.append(uKal)
    
    uddot_old = np.copy(uddot)
    uddot_old_true = uddot_true
    uddot = np.linalg.solve(m_eff,P_eff[0])
    uddot_true = np.linalg.solve(m_eff,P_eff_true[0])
    w1[0] = (np.random.rand()-.5)*.01
    u = u_pred+(dt**2)*Beta*uddot
    u_true = u_pred_true+(dt**2)*Beta*uddot_true
    udot = udot + dt*((1-gamma)*uddot_old+gamma*uddot)
    udot_true = udot_true + dt*((1-gamma)*uddot_old_true+gamma*uddot_true)

    u_pred = u+dt*udot+(dt**2)*(1-2*Beta)*uddot/2
    u_pred_true = u_true+dt*udot_true+(dt**2)*(1-2*Beta)*uddot_true/2
    #print(u_pred)
    U_hist.append(u)
    U_hist_true.append(u_true)
    
U_hist = np.vstack(U_hist)  
U_hist_true = np.vstack(U_hist_true)  
P_hist = np.vstack(P_hist)   
#U_Kalman = np.vstack(U_Kalman)

#plt.plot(U_hist[:,0]) 
plt.figure(21)
T = len(U_hist[:,0])*dt
time = np.arange(0,T,dt)

#plt.plot(time,U_Kalman[:,0],label = "Kalman Estimate") 
#plt.plot(time,U_hist[:,0],label = "Noisy Result") 
plt.plot(time,U_hist_true[:,0],'--',label = "True Result") 
plt.legend()

# import os
# import numpy as np
# import matplotlib.pyplot as plt

# calcDisp = os.path.join("/Users/juanmeriles/Google Drive (Backup)/Documents/UC Berkeley/Research"+"/Proj_LinTest3/", "IntegratorDisp.txt")
    
# analyticDisp=[]
# f=open(calcDisp)

# analyticDisp = []
# for line in f.readlines():
#     analyticDisp.append(float(line))

# f.close
# analyticDisp = analyticDisp[0:-2]
# plt.figure(2)
# plt.plot(analyticDisp,label = "Hybrid Result")
# plt.plot(U_hist.T[0],label = "Analytical Result")
# plt.legend()
#Dynamic analysis
#Model.Eig()
#Model.Create_C("Rayleigh",[1,1,1])
#An=Analyze()
#GM=GMread(Ground_Motion)
#[Model.Udy,Model.U_dotdy,Model.U_ddotdy]=An.Dyn_Newmark(Model,DynOp,[],GM)

#Static Analysis
#An.Displ(Model,el)


#plot=PlotModel()
#plot.Struct(Model,el)
#plot.DefShape(Model,el)
#animation = plot.AnimateDefShape(Model,el,100)
