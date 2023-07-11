#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 19:53:02 2020

@author: jim53
"""


import numpy as np
import sympy as sym
from Element import Element


class Analyze:
    
    def Displ(self,Model,el):
        Model.Create_Loading(el)
        Model.U=np.linalg.solve(Model.K,Model.P)

        Model.Ut=np.pad(Model.U,((0,len(Model.BDOF)),(0,0)),'constant')
        
        #Post processing-> turn the U->v->q->P check global equilibrium
        
        if Model.ModelSet[1] =='line':
            Model.v=np.matmul(np.transpose(Model.Trf),Model.U)
            Model.q=np.matmul(Model.Ks,Model.v)+np.transpose(Model.Q0L)
            Model.P_eq=(np.matmul(Model.Trt,Model.q))
            Model.P_r=(Model.P_eq-(Model.Par))
        
    def Dyn_Newmark(self,Model,DynOp,P,GM):
        #DynOp is a vector made up of [beta,gamma,dt,dir]
        #dir should be of the form [1,1,1] where 1 indicates the GM is applied in the x,y,z dir respectively
        #P can be an empty vector or it can be a dynamic load matrix (Will create)
        #a code to make simple loading functions like ramps
        #Should be (FDOF number of rows, Time/dt number of columns)
        #GM should be acceleration ground motion data in one row
        flag = 0

        if P==[]:
            P=np.zeros((len(Model.FDOF),GM.shape[0]))
            g=386.4 #in/s^2, this is the standard but I should probably let people change it
            GM=386.4*GM
            LoadingLen = GM.shape[0]
            flag = 1
        if GM==[]:
            GM=np.zeros((1,P.shape[1]))
            LoadingLen = P.shape[1]
            Res = np.zeros((1,P.shape[0]))
        
        
        beta=DynOp[0]
        gamma=DynOp[1]
        dt=DynOp[2]
        c0=1/(beta*dt**2)
        c1=1/(beta*dt)
        c2=gamma/(beta*dt)
        c3=(1/(2*beta)-1)
        c4=(gamma/beta)-1
        c5=dt*(gamma/(2*beta)-1)
        
        Keff=c0*Model.M+c2*Model.C+Model.K
        U=np.zeros((len(Model.FDOF),LoadingLen))
        U_dot=np.zeros((len(Model.FDOF),LoadingLen))
        U_ddot=np.zeros((len(Model.FDOF),LoadingLen))
        Pbar = np.zeros((len(Model.FDOF),LoadingLen))
        
        BigK=np.matmul(Model.Trt,np.matmul(Model.Ks,np.transpose(Model.Trt)))
        
        K=Model.K
        C=Model.C
        M=Model.M

        
        #Set of influence vectors for each fixed dof
        #R=-np.linalg.solve(BigK[0:len(Model.FDOF),0:len(Model.FDOF)],BigK[0:len(Model.FDOF),len(Model.FDOF):Model.NDOF])
        #print(R)
        #MR=np.matmul(M,R)
    
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
        # if Model.ModelSet[0]=='Planar':
        #     for i in range(len(Model.BDOF)):
        #         if Model.BDOF[i]+3 % 3 == 0 and DynOp[3][0]==1:
        #             influence.append(1)
        #         elif (Model.BDOF[i]+3-1) % 3==0 and DynOp[3][1]==1:
        #             influence.append(1)
        #         else:
        #             influence.append(0)
                    
        # else:
        #     for i in range(len(Model.BDOF)):
        #         if Model.BDOF[i]+6 % 6 == 0 and DynOp[3][0]==1:
        #             influence.append(1)
        #         elif (Model.BDOF[i]+6-1) % 6==0 and DynOp[3][1]==1:
        #             influence.append(1)
        #         elif (Model.BDOF[i]+6-2) % 6==0 and DynOp[3][2]==1:
        #             influence.append(1)
        #         else:
        #             influence.append(0)
        influence=np.matrix(influence)   
        
        print(influence)
        if flag == 1:
            #P0=P[:,0]-np.transpose(np.matmul(MR,GM[0,0]*np.transpose(influence)))
            P0=P[:,0]-np.transpose(np.matmul(M,GM[0]*np.transpose(influence)))
        
            U_ddot[:,0]=np.transpose(np.linalg.solve(M,np.transpose(P0[0]-np.matmul(C,(U_dot[:,0]))-np.matmul(K,U[:,0]))))

            
            T=dt*GM.shape[0]
            t=np.linspace(0,T,GM.shape[0]-1)
            

            for i in range(GM.shape[0]-1):
                #Pbar=P[:,i+1]-np.transpose(np.matmul(MR,GM[0,i+1]*np.transpose(influence)))
                Pbar[:,i+1]=P[:,i+1]-np.transpose(np.matmul(M,GM[i+1]*np.transpose(influence)))
                #print(M)

                Peff=Pbar[:,i+1]+\
                    np.matmul((c0*M+c2*C),U[:,i])+np.matmul((c1*M+c4*C),U_dot[:,i])+\
                    np.matmul((c3*M+c5*C),U_ddot[:,i])
                
                U[:,i+1]=np.transpose(np.linalg.solve(Keff,np.transpose(Peff)))
                
                
                U_dot[:,i+1]=c2*(U[:,i+1]-U[:,i])-c4*U_dot[:,i]-c5*U_ddot[:,i]
                U_ddot[:,i+1]=c0*(U[:,i+1]-U[:,i])-c1*U_dot[:,i]-c3*U_ddot[:,i]
        else:
            P0 = P[:,0]
            U_ddot[:,0]=np.transpose(np.linalg.solve(M,np.transpose(P0[0]-np.matmul(C,(U_dot[:,0]))-np.matmul(K,U[:,0]))))
            T=dt*P.shape[1]
            t=np.linspace(0,T,P.shape[1]-1)
            
            for i in range(P.shape[1]-1):
                Pbar=P[:,i+1].T
                #print(Pbar)
                
                
                Peff=Pbar+\
                    np.matmul((c0*M+c2*C),U[:,i])+np.matmul((c1*M+c4*C),U_dot[:,i])+\
                    np.matmul((c3*M+c5*C),U_ddot[:,i])
                
                print(Peff)
                U[:,i+1]=np.transpose(np.linalg.solve(Keff,np.transpose(Peff)))
                
                
                U_dot[:,i+1]=c2*(U[:,i+1]-U[:,i])-c4*U_dot[:,i]-c5*U_ddot[:,i]
                U_ddot[:,i+1]=c0*(U[:,i+1]-U[:,i])-c1*U_dot[:,i]-c3*U_ddot[:,i]
                
                Res = Keff@U[:,i+1]
                #print(Res)
        return[U,U_dot,U_ddot,Pbar]
    

    # def Hyb_alphaOS(self,Model,DynOp,P,GM):       
    #     u = np.zeros(len(Model.U))
    #     udot = np.zeros(len(Model.U))
    #     uddot = np.zeros(len(Model.U))
    #     u_pred = np.zeros(len(Model.U))
        
    #     m_eff = Model.M + (1-alpha)*dt*gamma*Model.C + dt**2*Beta*(1-alpha)*Model.K
    #     u_pred = u+dt*udot+(dt**2)*(1-2*Beta)*uddot/2
        
        
    #     P_resisting = Model.K@u_pred
        
    #     P_tot = Pm+P_resisting
        
    #     p_eff = (1-alpha)*P_next+alpha*P-(1-alpha)*P_tot-alpha*P_old- ((1-alpha)*C*dt*(1-gamma)+alpha*(dt**2)*Beta*K)*uddot-C*udot
    #     uddot_old = uddot
    #     uddot = np.linalg.solve(m_eff,p_eff)
    #     u = u_pred+(dt**2)*Beta*uddot
    #     udot = udot + dt*((1-gamma)*uddot_old+gamma*uddot)
    #     u_pred = u+dt*udot+(dt**2)*(1-2*Beta)*uddot/2
        
    #     return u,udot,uddot,u_pred,Model
        
        
        
        
                
        
        
        
        
        