#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:14:42 2020

@author: jim53
"""

import numpy as np
import sympy as sym
from Element import Element
import scipy.linalg


                 
class Model:
    
    def __init__(self, NODE, CON,Node_Mass,selfweight,ModelSet):
        self.NODE = NODE
        self.Node_Mass=Node_Mass
        self.CON = CON
        self.BOUN = []
        self.LOAD = []
        self.DOF = []
        self.numel = len(CON)
        self.selfweight=selfweight
        self.ModelSet=ModelSet
        self.hybridDOFS = []
        self.hybrid = 0
        if len(ModelSet) == 2:
            self.ModelSet.append("NULL")
        count=1
        self.DynOp = []
        #number of forces per element
        #maybe this should be a property of each element and not of the model

        
        if ModelSet[0]=='Planar' and ModelSet[1] == 'line':
            self.nf=6
        elif (ModelSet[0] == '3D' and ModelSet[1] == 'line'):
            self.nf=12
        elif (ModelSet[0] == 'Planar' and (ModelSet[1] == 'FEA_Tri')):
            self.nf = 6
        elif (ModelSet[0] == 'Planar' and ModelSet[1] == 'FEA_4nQuad'):
            self.nf = 8
        elif (ModelSet[0] == 'Planar' and ModelSet[1] == 'FEA_8nQuad'):
            self.nf = 16
        
        
        
        #Initializes the BOUNDS LOADS and generates global DOFS
        
        if ModelSet[0]=='Planar' and ModelSet[1] == 'line':
            for i in range(len(NODE)):
                self.LOAD.append([0,0,0])
                self.BOUN.insert(i,[0,0,0])
                self.DOF.insert(i,[count,count+1,count+2])
                count=count+3
            self.NDOF=count-1
        elif ModelSet[0] =='3D' and ModelSet[1] == 'line':
            for i in range(len(NODE)):
                self.LOAD.append([0,0,0,0,0,0])
                self.BOUN.insert(i,[0,0,0,0,0,0])
                self.DOF.insert(i,[count,count+1,count+2,count+3,count+4,count+5])
                count=count+6
            self.NDOF=count-1
        elif (ModelSet[0] =='Planar' ):
            for i in range(len(NODE)):
                self.LOAD.append([0,0])
                self.BOUN.insert(i,[0,0])
                self.DOF.insert(i,[count,count+1])
                count=count+2
            self.NDOF=count-1
        
                      
        #Initialize Elements
    def init_el(self):
        el=[]
        for i in range(self.numel):
            #check element type, do different for truss and frame
            
            if self.ModelSet[0]=='Planar' and self.ModelSet[1] =='line':
                el.insert(i,Element(['Bernoulli','Planar']))
                el[i].node=[self.NODE[self.CON[i][0]-1],self.NODE[self.CON[i][1]-1]]
                el[i].get_Length()
                el[i].BOUN=[self.BOUN[self.CON[i][0]-1],self.BOUN[self.CON[i][1]-1]]
            elif self.ModelSet[0] == '3D' and self.ModelSet[1] == 'line':
                el.insert(i,Element(['Bernoulli','3D']))
                el[i].node=[self.NODE[self.CON[i][0]-1],self.NODE[self.CON[i][1]-1]]
                el[i].get_Length()
                el[i].BOUN=[self.BOUN[self.CON[i][0]-1],self.BOUN[self.CON[i][1]-1]]
            elif self.ModelSet[0] =='Planar' and self.ModelSet[1] == 'FEA_Tri':
                el.insert(i,Element(['Tri','Planar']))
                el[i].node=[self.NODE[self.CON[i][0]-1],self.NODE[self.CON[i][1]-1],self.NODE[self.CON[i][2]-1]]
                el[i].get_Area()
                el[i].BOUN=[self.BOUN[self.CON[i][0]-1],self.BOUN[self.CON[i][1]-1],self.BOUN[self.CON[i][2]-1]]
                el[i].w = np.array([0,0,0,0,0,0])
            elif self.ModelSet[0] =='Planar' and self.ModelSet[1] == 'FEA_4nQuad':
                el.insert(i,Element(['4nQuad','Planar']))
                el[i].node=[self.NODE[self.CON[i][0]-1],self.NODE[self.CON[i][1]-1],self.NODE[self.CON[i][2]-1],self.NODE[self.CON[i][3]-1]]
                el[i].get_Area()
                el[i].BOUN=[self.BOUN[self.CON[i][0]-1],self.BOUN[self.CON[i][1]-1],self.BOUN[self.CON[i][2]-1],self.BOUN[self.CON[i][3]-1]]
                el[i].w = np.array([0,0,0,0,0,0,0,0])
            elif self.ModelSet[0] =='Planar' and self.ModelSet[1] == 'FEA_8nQuad':
                el.insert(i,Element(['8nQuad','Planar']))
                el[i].node=[self.NODE[self.CON[i][0]-1],self.NODE[self.CON[i][1]-1],self.NODE[self.CON[i][2]-1],self.NODE[self.CON[i][3]-1],\
                            self.NODE[self.CON[i][4]-1],self.NODE[self.CON[i][5]-1],self.NODE[self.CON[i][6]-1],self.NODE[self.CON[i][7]-1]]
                el[i].get_Area()
                el[i].BOUN=[self.BOUN[self.CON[i][0]-1],self.BOUN[self.CON[i][1]-1],self.BOUN[self.CON[i][2]-1],self.BOUN[self.CON[i][3]-1],\
                            self.BOUN[self.CON[i][4]-1],self.BOUN[self.CON[i][5]-1],self.BOUN[self.CON[i][6]-1],self.BOUN[self.CON[i][7]-1]]
                
            
            el[i].CON=self.CON[i]
            
            if self.ModelSet[0]=='Planar' and self.ModelSet[1] =='line':
                el[i].REL=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
                el[i].J=0
            elif self.ModelSet[0] == '3D' and self.ModelSet[1] =='line':
                el[i].REL=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
            if self.ModelSet[1] == 'line':
                el[i].q0=0
                el[i].t0=0
                el[i].w_x=0
                el[i].w_z=0
                el[i].w_y=0
                el[i].m_x=0
                el[i].m_y=0
                el[i].m_z=0
                el[i].rho=0
                el[i].hybrid = 0
            
        return el
    
    def elMat(self,el):
        #creates the element Stiffness and Flexibility matrices
        gl=0
        self.Q0L=[]
        if self.ModelSet[1] == 'line':
            for i in range (self.numel):
                 #Updates Boundaries of the element
                 A=el[i].A
                 L=el[i].L
                 E=el[i].E
                 Iz=el[i].Iz
                 
                 if self.ModelSet[0]!='Planar':
                     G=el[i].G
                     J=el[i].J
                     Iy=el[i].Iy
                 
                 if self.selfweight==1:
                    [el[i].ROT,br]=el[i].Create_Rot()
                    #swg is selfweight global
                    el[i].swg=np.matrix([0,el[i].rho,0])  
                    #swl is selfweight local
                    el[i].swl=np.linalg.solve(el[i].ROT,np.transpose(el[i].swg))
                    el[i].w_x=el[i].swl[0,0]+el[i].w_x
                    el[i].w_y=el[i].swl[1,0]+el[i].w_y
                    if self.ModelSet[0]!='Planar':
                        el[i].w_z=el[i].swl[2,0]+el[i].w_z
                
                 if self.ModelSet[0]=='Planar':
                     
                     if el[i].hybrid == 0 and el[i].custom == 0:                   
                         el[i].ke=np.matrix([[E*A/L,0,0],\
                                             [0,4*E*Iz/L,2*E*Iz/L],\
                                             [0,2*E*Iz/L,4*E*Iz/L]])
                         el[i].fe=np.matrix([[L/(E*A),0,0],\
                                             [0,L/(3*E*Iz),-L/(6*E*Iz)],\
                                             [0,-L/(6*E*Iz),L/(3*E*Iz)]])
                     else:
                         #print(el[i].ke)
                     #else:            
                         #boun = np.hstack(el[i].BOUN)
                         # if (el[i].controlDOF[0] == 0 and el[i].controlDOF[3] == 0) and (boun[0] == 1 or boun[3] == 1):
                         #    el[i].ke=np.matrix([[10000000000,0,0],\
                         #                         [0,0,0],\
                         #                         [0,0,0]])
                         #    el[i].fe=np.matrix([[0,0,0],\
                         #                         [0,0,0],\
                         #                         [0,0,0]])
                         # else:
                         #    el[i].ke=np.matrix([[0,0,0],\
                         #                         [0,0,0],\
                         #                         [0,0,0]])
                         #    el[i].fe=np.matrix([[0,0,0],\
                         #                         [0,0,0],\
                         #                         [0,0,0]])
                         self.hybrid = 1
                         
                                     
                         eldofs = np.hstack([self.DOF[el[i].CON[0]-1],self.DOF[el[i].CON[1]-1]])
                         self.hybridDOFS.append(eldofs)
                         self.interfaceDOFS = []
                         for p in range(len(eldofs)):
                             for o in range(len(el)):
                                 if el[o].hybrid == 0:
                                     othereldofs = np.hstack([self.DOF[el[o].CON[0]-1],self.DOF[el[o].CON[1]-1]])   
                                     if eldofs[p] in othereldofs and eldofs[p] not in self.interfaceDOFS:
                                        self.interfaceDOFS.append(eldofs[p]-1)
                                        
                        
                    
                     Jmat=[el[i].Iz]
                     el[i].create_Shape()
                     el[i].me=el[i].rho*(el[i].A*np.matmul(np.transpose(el[i].Nu),\
                            el[i].Nu)+np.matmul(np.transpose(el[i].Nt),\
                            np.matmul(Jmat,el[i].Nt)))
                 else:
                     if el[i].hybrid == 0 and el[i].custom == 0:
                         el[i].ke=np.matrix([[E*A/L,0,0,0,0,0],\
                                             [0,G*J/L,0,0,0,0],\
                                             [0,0,4*E*Iy/L,2*E*Iy/L,0,0],\
                                             [0,0,2*E*Iy/L,4*E*Iy/L,0,0],\
                                             [0,0,0,0,4*E*Iz/L,2*E*Iz/L],\
                                             [0,0,0,0,2*E*Iz/L,4*E*Iz/L]])
                         el[i].fe=np.matrix([[L/(E*A),0,0,0,0,0],\
                                             [0,L/(G*J),0,0,0,0],\
                                             [0,0,L/(3*E*Iy),-L/(6*E*Iy),0,0],\
                                             [0,0,-L/(6*E*Iy),L/(3*E*Iy),0,0],\
                                             [0,0,0,0,L/(3*E*Iz),-L/(6*E*Iz)],\
                                             [0,0,0,0,-L/(6*E*Iz),L/(3*E*Iz)]])
                     else:
                         el[i].ke=np.matrix([[0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0]])
                         el[i].fe=np.matrix([[0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0],\
                                             [0,0,0,0,0,0]])
                         self.hybrid = 1
                         
                         
                         
                         eldofs = np.hstack([self.DOF[el[i].CON[0]-1],self.DOF[el[i].CON[1]-1]])
                         self.hybridDOFS.append(eldofs)
                         self.interfaceDOFS = []
                         for p in range(len(eldofs)):
                             for o in range(len(el)):
                                 if el[o].hybrid == 0:
                                     othereldofs = np.hstack([self.DOF[el[o].CON[0]-1],self.DOF[el[o].CON[1]-1]])       
                                     if eldofs[p] in othereldofs and eldofs[p] not in self.interfaceDOFS:
                                        self.interfaceDOFS.append(eldofs[p]-1)
                                             
                         
                         
                         
                     Jmat=np.matrix([[el[i].J,0,0],\
                                  [0,el[i].Iy,0],\
                                  [0,0,el[i].Iz]])
                     
                 
                     el[i].create_Shape()
                     
                     # el[i].me=el[i].rho*(el[i].A*np.matmul(np.transpose(el[i].Nu),\
                     #        el[i].Nu)+np.matmul(np.transpose(el[i].Nt),\
                     #        np.matmul(Jmat,el[i].Nt)))
    ### Important
                     #el[i].me = el[i].rho*((el[i].A*el[i].Nu.T@el[i].Nu)+el[i].Nt.T@Jmat@el[i].Nt)
                     el[i].me = el[i].rho*((el[i].A*el[i].Nu.T@el[i].Nu))
                    
                 
                 
                 
                 #ep is my local coordinate that just goes to the length of the element
                 ep=np.linspace(0,el[i].L,1000)
                 x=sym.symbols('x')
                 for k in range(self.nf):
                     for j in range(self.nf):
                         if el[i].me[k,j]==0:
                             el[i].me[k,j]=0.0
                         else:
                             f=sym.lambdify(x,el[i].me[k,j],"numpy")
                             el[i].me[k,j]=np.trapz(f(ep),ep)
                 el[i].me=el[i].me.astype('float64')
                 
                 
                 
                 
                 el[i].create_LVec()
                 el[i].release()
                 
                 el[i].create_Q0()
                 #el[i].create_v0()
                
                 if self.ModelSet[0]=='Planar':
                    b=np.matrix([[-1,0,0],\
                                  [0,1/L,1/L],\
                                  [0,1,0],\
                                  [1,0,0],\
                                  [0,-1/L,-1/L],\
                                  [0,0,1]])
                 else: 
                    b=np.matrix([[-1,0,0,0,0,0],\
                     [0,0,0,0,1/el[i].L,1/el[i].L],\
                     [0,0,-1/el[i].L,-1/el[i].L,0,0],\
                     [0,-1,0,0,0,0],\
                     [0,0,1,0,0,0],\
                     [0,0,0,0,1,0],\
                     [1,0,0,0,0,0],\
                     [0,0,0,0,-1/el[i].L,-1/el[i].L],\
                     [0,0,1/el[i].L,1/el[i].L,0,0],\
                     [0,1,0,0,0,0],\
                     [0,0,0,1,0,0],\
                     [0,0,0,0,0,1]])
                        
                 self.Q0L.append(np.transpose(np.matmul(b,np.transpose(el[i].eQ0))))
            self.Q0L=np.hstack(self.Q0L) 
            
            if self.hybridDOFS != []:
                self.hybridDOFS = np.hstack(self.hybridDOFS)-1 
                self.Ks_hybrid=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Ks=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Ms=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Fs=np.zeros((self.nf*self.numel,self.nf*self.numel))
            
            
            #create diagonal total stiffness/flex
            g=0
            for n in range(self.numel):
                for j in range (len(el[n].ke)):
                    for i in range(len(el[n].ke)):
                        self.Ks[g+i][g+j]=el[n].ke[i,j]
                        self.Ms[g+i][g+j]=el[n].me[i,j]
                        if self.hybrid == 1:
                            if el[n].hybrid == 1:
                                self.Ks_hybrid[g+i][g+j] = 0
                            else:
                                self.Ks_hybrid[g+i][g+j] = el[n].ke[i,j]
                        # Fs[g+i][g+j]=el[n].fe[i][j]
                g=g+j+1;
                         
        elif self.ModelSet[1] == 'FEA_Tri':
            
            self.Ks=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Ms=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Fs=np.zeros((self.nf*self.numel,self.nf*self.numel))
            
            
            if self.ModelSet[2] =='Stress':
                for i in range (self.numel):
                    el[i].create_Shape()
                    el[i].D = (el[i].E)/(1-el[i].v**2)*np.matrix([[1,el[i].v,0],\
                                                        [el[i].v,1,0],\
                                                        [0,0,(1-el[i].v)/2]])
                    el[i].ke = el[i].th*el[i].Be.T @ el[i].D @ el[i].Be
                
                                   
                        
            elif self.ModelSet[2] =='Strain':
                for i in range (self.numel):      
                    el[i].create_Shape()
                    el[i].D = el[i].E*(1-el[i].v)/((1+el[i].v)*(1-2*el[i].v))*np.matrix([[1,el[i].v/(1-el[i].v),0],\
                                                                               [el[i].v/(1-el[i].v),1,0],\
                                                                               [0,0,(1-2*el[i].v)/(2*(1-el[i].v))]])
                    el[i].ke = el[i].th*el[i].Be.T @ el[i].D @ el[i].Be
                    
            g=0
            for n in range(self.numel):
                for j in range (len(el[n].ke)):
                    for i in range(len(el[n].ke)):
                        self.Ks[g+i][g+j]=el[n].ke[i,j]
                        #self.Ms[g+i][g+j]=el[n].me[i,j]
                        # Fs[g+i][g+j]=el[n].fe[i][j]
                g=g+j+1;
                
        elif self.ModelSet[1] == 'FEA_4nQuad':
            
            self.Ks=np.zeros((self.nf*self.numel,self.nf*self.numel))
            self.Ms=np.zeros((self.nf*self.numel,self.nf*self.numel))
            #self.Fs=np.zeros((self.nf*self.numel,self.nf*self.numel))
            ek=[-0.57735,0.57735]
            nl=[-0.57735,0.57735]
            
            e = sym.symbols('e')
            n = sym.symbols('n')
            
            if self.ModelSet[2] =='Stress':
                for i in range (self.numel):
                    el[i].ke = np.zeros((self.nf,self.nf))
                    el[i].create_Shape()
                    el[i].D = (el[i].E)/(1-el[i].v**2)*np.array([[1,el[i].v,0],\
                                                        [el[i].v,1,0],\
                                                        [0,0,(1-el[i].v)/2]])
                        
                    #use gauss quadrature to get stiffness matrix
                    
                    #el[i].ke = el[i].th*el[i].Be.T @ el[i].D @ el[i].Be
                    for k in range(len(ek)): 
                        for l in range(len(nl)):
                            temp = (el[i].Be.T @ el[i].D @ el[i].Be)
                            for g in range(self.nf):
                                for h in range(self.nf):

                                    temp[g][h] = float(temp[g][h].subs(e,ek[k]).subs(n,nl[l]))
                                    
                            el[i].ke = el[i].ke+temp
                
                                   
                        
            elif self.ModelSet[2] =='Strain':
                for i in range (self.numel):    
                    el[i].ke = np.zeros((self.nf,self.nf))
                    el[i].create_Shape()
                    el[i].D = el[i].E*(1-el[i].v)/((1+el[i].v)*(1-2*el[i].v))*np.array([[1,el[i].v/(1-el[i].v),0],\
                                                                               [el[i].v/(1-el[i].v),1,0],\
                                                                               [0,0,(1-2*el[i].v)/(2*(1-el[i].v))]])
                    for k in range(len(ek)): 
                        for l in range(len(nl)):
                            temp = (el[i].Be.T @ el[i].D @ el[i].Be)
                            for g in range(self.nf):
                                for h in range(self.nf):
                                    temp[g][h] = float(temp[g][h].subs(e,ek[k]).subs(n,nl[l]))

                            el[i].ke = el[i].ke+temp

                            
                        
                    
                    
            g=0
            for n in range(self.numel):
                for j in range (len(el[n].ke)):
                    for i in range(len(el[n].ke)):
                        self.Ks[g+i][g+j]=el[n].ke[i,j]
                        #self.Ms[g+i][g+j]=el[n].me[i,j]
                        # Fs[g+i][g+j]=el[n].fe[i][j]
                g=g+j+1;
        return el
        

        
    def Create_B(self,el):
        if self.ModelSet[1] =='line':
            Cindex=0 #track the Global index for which forces are Constrained
            CForces=[]
            
            bg=np.zeros(self.numel)
            bg=bg.tolist()
            Ab=np.zeros(self.numel)
            Ab=Ab.tolist()
            self.br=[]
            for i in range(self.numel):
                 #Creates b matrix for each element that changes the 12 element forces
                 #into 6 basic element forces
                L=el[i].L
                if self.ModelSet[0]=='Planar':
                    b=np.matrix([[-1,0,0],\
                                  [0,1/L,1/L],\
                                  [0,1,0],\
                                  [1,0,0],\
                                  [0,-1/L,-1/L],\
                                  [0,0,1]])
                else: 
                    b=np.matrix([[-1,0,0,0,0,0],\
                     [0,0,0,0,1/el[i].L,1/el[i].L],\
                     [0,0,-1/el[i].L,-1/el[i].L,0,0],\
                     [0,-1,0,0,0,0],\
                     [0,0,1,0,0,0],\
                     [0,0,0,0,1,0],\
                     [1,0,0,0,0,0],\
                     [0,0,0,0,-1/el[i].L,-1/el[i].L],\
                     [0,0,1/el[i].L,1/el[i].L,0,0],\
                     [0,1,0,0,0,0],\
                     [0,0,0,1,0,0],\
                     [0,0,0,0,0,1]])
                
                [el[i].ROT,br]=el[i].Create_Rot()
                self.br.append(br)
                    
               
            
            #     creates matrices necessary for the Bf matrix
                bg[i]=np.matmul(br,b)
                el[i].id_v=self.DOF[el[i].CON[0]-1]+self.DOF[el[i].CON[1]-1]
                Ab[i]=np.zeros((self.nf,self.NDOF))
                
                
                
            #     This searches for forces that are unreleased "constrained" in each
            #     element. The released forces are removed from the Bf matrix
                for w in range(int(self.nf/2)):
                    if el[i].REL[w]!=0:
                        CForces.append(Cindex)  
                        Cindex=Cindex+1
                    else:
                        Cindex=Cindex+1
                
            #     This creates matrix Ab (Fillipou pg. 59, which is a mapping of which
            #     Degrees of freedom are involved in an element on each end. It is a 2d
            #     Version of the id matrix, showing how each local dof corresponds to a
            #     global dof
                
                for j in range(len(el[i].id_v)):
                    Ab[i][j][el[i].id_v[j]-1]=1
            
            # %This creates the Bf matrix
            #This creates the Transformation matrix Tr (As opposed to Bf, this matrix is if we use full element stiffness matrix)
            B=list()
            Tr=list()
            for i in range(self.numel):
                B.append(np.matmul(np.transpose(Ab[i]),bg[i]))
                Tr.append(np.matmul(np.transpose(Ab[i]),self.br[i]))
                
            
            B=np.hstack(B)
            Tr=np.hstack(Tr)
    
            # This removes the released columns from the bf matrix as the forces
            # associated with them is 0
            
            #B=B[:,CForces]
            #Why is this commented out, need to look into this
            
            
            #Rearrange Bf matrix to put boundary rows last
            count=0
            self.BDOF=[]
            self.FDOF=[]
            for n in range(len(self.BOUN)):
                    for j in range(int(self.nf/2)):
                        if self.BOUN[n][j]==1:
                            self.BDOF.append(count)                
                        else:
                            self.FDOF.append(count)
                        
                        count=count+1;
                

            freehdofs = []
            compdofs = []
            
            for i in range(len(self.hybridDOFS)):
                for j in range(len(self.FDOF)):
                    if self.FDOF[j] == self.hybridDOFS[i]:
                        freehdofs.append(self.FDOF[j])
                    
            for i in range(len(self.FDOF)):
                flag = 0
                for j in range(len(self.hybridDOFS)):
                    if self.FDOF[i] == self.hybridDOFS[j]:
                        flag = 1
                if flag == 0:
                    compdofs.append(self.FDOF[i])
            

            self.FDOF = freehdofs+compdofs
            
            
            # This is the "Total" Bf and transformation matrix, which includes the boundary DOFS
            self.Bt=np.vstack([B[self.FDOF,:],B[self.BDOF,:]])
            self.Trt=np.vstack([Tr[self.FDOF,:],Tr[self.BDOF,:]])
            
            # This only includes the free dofs, boundaries and releases removed
            self.Bf=B[self.FDOF,:]; 
            self.Trf=Tr[self.FDOF,:] #releases not removed for Assembly transformation matrix
            tempdoforder = np.hstack([self.FDOF,self.BDOF]).tolist()
            if self.hybrid == 1:
                for i in range(len(self.hybridDOFS)):
                    self.hybridDOFS[i] = tempdoforder.index(self.hybridDOFS[i])
                for i in range(len(self.interfaceDOFS)):
                    self.interfaceDOFS[i] = tempdoforder.index(self.interfaceDOFS[i])
                                

            #check this
        else:
            Cindex=0 #track the Global index for which forces are Constrained
            CForces=[]
            Ab=np.zeros(self.numel)
            Ab=Ab.tolist()
            for i in range(self.numel):
                if self.ModelSet[1] == 'FEA_Tri':
                    el[i].id_v=self.DOF[el[i].CON[0]-1]+self.DOF[el[i].CON[1]-1]+self.DOF[el[i].CON[2]-1]
                    r=self.nf/3
                elif self.ModelSet[1] == 'FEA_4nQuad':
                    el[i].id_v=self.DOF[el[i].CON[0]-1]+self.DOF[el[i].CON[1]-1]+self.DOF[el[i].CON[2]-1]+self.DOF[el[i].CON[3]-1]
                    r = self.nf/4
                Ab[i]=np.zeros((self.nf,self.NDOF))
                for j in range(len(el[i].id_v)):
                    Ab[i][j][el[i].id_v[j]-1]=1
                        
            B=list()
            Tr=list()
            for i in range(self.numel):
                B.append(np.transpose(Ab[i]))
                Tr.append((np.transpose(Ab[i])))
                
            
            B=np.hstack(B)
            Tr=np.hstack(Tr)
            
            count=0
            self.BDOF=[]
            self.FDOF=[]
            
            for n in range(len(self.BOUN)):
                    for j in range(int(r)):
                        if self.BOUN[n][j]==1:
                            self.BDOF.append(count)                
                        else:
                            self.FDOF.append(count)
                        
                        count=count+1;
                    
            self.test = Tr
            
            # This is the "Total" Bf and transformation matrix, which includes the boundary DOFS
            #Throws the Bound DOFS to the end of the matrix
            self.Bt=np.vstack([B[self.FDOF,:],B[self.BDOF,:]])
            self.Trt=np.vstack([Tr[self.FDOF,:],Tr[self.BDOF,:]])
            
            # This only includes the free dofs, boundaries and releases removed
            self.Bf=B[self.FDOF,:]; 
            self.Trf=Tr[self.FDOF,:] #releases not removed for Assembly transformation matrix
                     

                    
                        
            # this is the global stiffness matrix
    def Create_K(self):
        self.K=np.matmul(self.Trf,np.matmul(self.Ks,np.transpose(self.Trf)))
        
        if self.hybrid ==1:
            self.Kh=np.matmul(self.Trf,np.matmul(self.Ks_hybrid,np.transpose(self.Trf)))

        # for i in range(self.K.shape[0]):
        #     for j in range(self.K.shape[1]):
        #         if self.K[i,j] == 0: 
        #             self.K[i,j] = 10**-18
    
    def Create_M(self):
        self.M= self.Trf @ self.Ms @ self.Trf.T
        #np.matmul(self.Trf,np.matmul(self.Ms,np.transpose(self.Trf)))
        
        
        # we need to have some kind of rotational mass component in the diag
        # of the mass matrix. I am not yet sure what the rotational mass of 
        # a point mass is
        if any(self.Node_Mass!=np.zeros((len(self.NODE)))):
            self.Mn=np.zeros((np.size(self.DOF),np.size(self.DOF)))
            count=0
            for i in range(len(self.DOF)):
                for j in range(len(self.DOF[0])):
                    if (self.DOF[i][j]+self.nf/2-1) % (self.nf/2)==0:
                        #Mn stands for Mass node
                        self.Mn[count,count]=self.Node_Mass[i]
                    if (self.DOF[i][j]+self.nf/2-1-1) % (self.nf/2)==0:
                       self. Mn[count,count]=self.Node_Mass[i]
                    
                    if self.ModelSet[0]!='Planar':
                        
                        if (self.DOF[i][j]+self.nf/2-2-1) % (self.nf/2)==0:
                            self.Mn[count,count]=self.Node_Mass[i]
                            
                    
                    count=count+1
                    
            #Mnf stands for Mass node free
            self.Mnf=np.matrix([[self.Mn[i,j] for i in self.FDOF] for j in self.FDOF])
            self.M=self.M+self.Mnf
                
        for i in range(self.M.shape[0]):
            for j in range(self.M.shape[1]):
                if self.M[i,j] == 0:
                    self.M[i,j] = 10**-18
                    pass
    def Create_Loading(self,el):
        self.Pw=[]
        self.Pf=[]
        self.Pwf=[]
        
        if self.ModelSet[1] == 'line':
            for i in range(self.numel):
                
                #This is the element distributed loads converted to applied nodal loads in local coords
                #once this vector is created,it is rotated to global and the dofs are placed in the right
                #location with the Tr matrix
                if self.ModelSet[0]=='Planar':
                    Pwfe=[-el[i].w_x*el[i].L/2,-el[i].w_y*el[i].L/2,0,\
                      -el[i].w_x*el[i].L/2,-el[i].w_y*el[i].L/2,0]
                else:
                
                    Pwfe=[-el[i].w_x*el[i].L/2,-el[i].w_y*el[i].L/2,-el[i].w_z*el[i].L/2,0,0,0,\
                          -el[i].w_x*el[i].L/2,-el[i].w_y*el[i].L/2,-el[i].w_z*el[i].L/2,0,0,0]
                        
                self.Pwf.append(Pwfe)
            self.Pwf=np.hstack(self.Pwf)
            self.Pwf=np.matmul(self.Trt,self.Pwf)    
            
            for i in range(len(self.NODE)):
                for j in range(len(self.BOUN[0])):
                    self.Pf.append(self.LOAD[i][j])
                        
            
            self.Pf=np.matrix(self.Pf)
            self.Pf=np.concatenate((self.Pf[0,self.FDOF],self.Pf[0,self.BDOF]),axis=1)
            
            
            self.Par=np.transpose(self.Pf-self.Pwf)
            
            #This can maybe be deleted
            #Par=np.transpose(self.Pf-self.Pwf)
            
            
            self.Pf=np.transpose(self.Pf)
            
            if self.selfweight==1:
                #P_sw is force due to nodal self weight. The distributed is treated as distributed load in Model.elMat
                self.P_sw=np.zeros((len(self.Pf)))
                for j in range(len(self.Node_Mass)):
                    DOForder=np.hstack([self.FDOF,self.BDOF])
                    if self.Node_Mass[j]!=0:
                        self.P_sw[list(DOForder).index(self.DOF[j][1]-1)]=-self.Node_Mass[j]*386.4
                    #Replace 386.4 with g eventually
                self.Pf=self.Pf+np.transpose(np.matrix(self.P_sw))
                
                
            
            for i in range(self.numel):  
                self.Pw.append(el[i].p_w)
            #Check on this
            self.Pw=np.hstack(self.Pw)
            self.Pw=np.matrix(self.Pw)
            self.Pw=np.matmul(self.Trt,np.transpose(self.Pw))
            self.Pt=self.Pf+self.Pw
            self.P=self.Pt[range(len(self.FDOF))]
            
        if self.ModelSet[1] == 'FEA_Tri':
            self.P_sw=[]
            for i in range(self.numel):
                #do distributed loading here, surface tractions
                count = 0
                Pwfe = np.array([[0,0,0,0,0,0]])
                for j in range(3):
                    t = np.array([el[i].w[count],el[i].w[count+1]])
                    
                    if j == 0:
                        L = np.sqrt((el[i].node[1][0]-el[i].node[0][0])**2+(el[i].node[1][1]-el[i].node[0][1])**2)
                        x = np.linspace(-1,1,1000)
                        N1 = (1/2)*(1+x)
                        N2 = (1/2)*(1-x)
                        N3 = 0
                        N = np.array([[np.trapz(N1,x),0,np.trapz(N2,x),0,0,0],\
                             [0,np.trapz(N1,x),0,np.trapz(N2,x),0,0]])
                        
                        Pwfe = Pwfe+(L/2)*N.T @ t

                        
                    if j == 1:
                        L = np.sqrt((el[i].node[2][0]-el[i].node[1][0])**2+(el[i].node[2][1]-el[i].node[1][1])**2)
                        x = np.linspace(-1,1,1000)
                        N1 = 0
                        N2 = (1/2)*(1+x)
                        N3 = (1/2)*(1-x)
                        
                        N = np.array([[0,0,np.trapz(N2,x),0,np.trapz(N3,x),0],\
                                      [0,0,0,np.trapz(N2,x),0,np.trapz(N3,x)]])
                        
                        
                        Pwfe = Pwfe+(L/2)*N.T @ t

                        
                    if j == 2:
                        L = np.sqrt((el[i].node[0][0]-el[i].node[2][0])**2+(el[i].node[0][1]-el[i].node[2][1])**2)
                        x = np.linspace(-1,1,1000)
                        N1 = (1/2)*(1-x)
                        N2 = 0
                        N3 = (1/2)*(1+x)
                        
                        N = np.array([[np.trapz(N1,x),0,0,0,np.trapz(N3,x),0],\
                                      [0,np.trapz(N1,x),0,0,0,np.trapz(N3,x)]])
                        
                        Pwfe = Pwfe+(L/2)*N.T @ t
                        
                    
                    count = count+2
                    #print(Pwfe)
                self.Pwf.append(Pwfe) 
                
            self.Pwf=np.hstack(self.Pwf)
            #Need to rearrange order here
            self.Pwf = self.Trt@self.Pwf[0]
            
                    
            #This is just point loads

            for i in range(len(self.NODE)):
                for j in range(len(self.BOUN[0])):
                    self.Pf.append(self.LOAD[i][j])
                    
            self.Pf=np.matrix(self.Pf)
            self.Pf=np.concatenate((self.Pf[0,self.FDOF],self.Pf[0,self.BDOF]),axis=1)
            self.Pt = self.Pf+self.Pwf
            #Here we put the selfweight

            if self.selfweight==1:
                #This is a body force
                for i in range(self.numel):
                    sw = np.array([0,-el[i].rho*386.4])
                    self.P_sw.append(el[i].Nu_int.T@sw)
                self.P_sw=np.hstack(self.P_sw)
                self.P_sw = self.Trt@self.P_sw
                self.Pt = self.Pt+self.P_sw
            self.P = self.Pt[0,range(len(self.FDOF))]
            self.P=(self.P.T)
            
        if self.ModelSet[1] == 'FEA_4nQuad':
            ek=[-0.57735,0.57735]
            nl=[-0.57735,0.57735]
            
            e = sym.symbols('e')
            n = sym.symbols('n')
            self.P_sw=[]

            for i in range(self.numel):
                #do distributed loading here, surface tractions
                count = 0
                x = np.linspace(-1,1,1000)
                
                Pwfe = np.array([[0,0,0,0,0,0,0,0]])
                for j in range(4):
                    t = np.array([el[i].w[count],el[i].w[count+1]])
                    
                    
                    if j == 0: #n is -1        
                        j1=float(el[i].x[0][0].diff(e).subs(n,-1))
                        j2=float(el[i].x[1][0].diff(e).subs(n,-1))
                        J = np.sqrt(j1**2+j2**2)
                        N1 = (1/2)*(1-x)
                        N2 = (1/2)*(1+x)
                        N3 = 0
                        N4 = 0
                        
                        N = np.array([[np.trapz(N1,x),0,np.trapz(N2,x),0,0,0,0,0],\
                                      [0,np.trapz(N1,x),0,np.trapz(N2,x),0,0,0,0]])
                                                
                    if j == 1: # e is 1
                        j1=float(el[i].x[0][0].diff(n).subs(e,1))
                        j2=float(el[i].x[1][0].diff(n).subs(e,1))
                        J = np.sqrt(j1**2+j2**2)
                        N1 = 0
                        N2 = (1/2)*(1-x)
                        N3 = (1/2)*(1+x)
                        N4 = 0
                        
                        
                        N = np.array([[0,0,np.trapz(N2,x),0,np.trapz(N3,x),0,0,0],\
                                      [0,0,0,np.trapz(N2,x),0,np.trapz(N3,x),0,0]])                        
                        
                        
                    if j == 2: #n is 1
                        j1=float(el[i].x[0][0].diff(e).subs(n,1))
                        j2=float(el[i].x[1][0].diff(e).subs(n,1))
                        J = np.sqrt(j1**2+j2**2)
                        N1 = 0
                        N2 = 0
                        N3 = (1/2)*(1+x)
                        N4 = (1/2)*(1-x)
                         
                        N = np.array([[0,0,0,0,np.trapz(N3,x),0,np.trapz(N4,x),0],\
                                      [0,0,0,0,0,np.trapz(N3,x),0,np.trapz(N4,x)]])                       
                        
                    if j == 3:# e is -1
                        j1=float(el[i].x[0][0].diff(n).subs(e,-1))
                        j2=float(el[i].x[1][0].diff(n).subs(e,-1))
                        J = np.sqrt(j1**2+j2**2)
                        N1 = (1/2)*(1-x)
                        N2 = 0
                        N3 = 0
                        N4 = (1/2)*(1+x)
                        
                        
                        N = np.array([[np.trapz(N1,x),0,0,0,0,0,np.trapz(N4,x),0],\
                                      [0,np.trapz(N1,x),0,0,0,0,0,np.trapz(N4,x)]])

                        
                    Pwfe = Pwfe+(J)*N.T @ t
                        
                    
                    count = count+2
                    #print(Pwfe)
                self.Pwf.append(Pwfe) 
                
            self.Pwf=np.hstack(self.Pwf)
            #Need to rearrange order here
            self.Pwf = self.Trt@self.Pwf[0]
            
                    
            #This is just point loads

            for i in range(len(self.NODE)):
                for j in range(len(self.BOUN[0])):
                    self.Pf.append(self.LOAD[i][j])
                    
            self.Pf=np.matrix(self.Pf)
            self.Pf=np.concatenate((self.Pf[0,self.FDOF],self.Pf[0,self.BDOF]),axis=1)
            self.Pt = self.Pf+self.Pwf
            #Here we put the selfweight
       
            if self.selfweight==1:
                #This is a body force
                for i in range(self.numel):
                    sw = np.array([0,-el[i].rho*386.4])
                    Jdet = abs(el[i].J[0][0]*el[i].J[1][1]-(el[i].J[1][0]*el[i].J[0][1]))
                    #unsure about the abs
                    tempsym = el[i].Nu.T @ sw * Jdet
                    temp = np.zeros((len(tempsym)))
                    for  j in range(len(temp)):
                        for k in range(len(ek)):
                            for l in range(len(nl)):
                                if tempsym[j]!=0:
                                    temp[j] = temp[j]+tempsym[j].subs(e,ek[k]).subs(n,nl[l])
                                   
                    self.P_sw.append(temp)
                self.P_sw=np.hstack(self.P_sw)
                self.P_sw = self.Trt@self.P_sw
                self.Pt = self.Pt+self.P_sw
            self.P = self.Pt[0,range(len(self.FDOF))]
            self.P=(self.P.T)
            
    
    def Eig(self):
        [self.w,self.vr]=scipy.linalg.eig(self.K,self.M)
        self.w=np.sqrt(self.w)
    def Create_C(self,DampType,DampOp):
        #DampType should be either Rayleigh or Caughy and Damp op should be
        #[percent(in percent eg. 5%=5),mode1,mode2]] for Rayleigh
        if DampType=="Rayleigh":
            DampOp[0] = DampOp[0]
            a0=np.real((DampOp[0]*2*self.w[DampOp[1]]*self.w[DampOp[2]])\
                /(self.w[DampOp[1]]+self.w[DampOp[2]]))
            print(a0)
            a1=np.real(2*DampOp[0]/(self.w[DampOp[1]]+self.w[DampOp[2]]))
            print(a1)
            self.C=a0*self.M+a1*self.K
            return [a0,a1]
        # if DampType=='Caughey':
            
    def CheckHybridEl(self,el):
        #This should work for what we're doing in 2D but 3D is gonna be a bit tougher
        #need a better solution
        print("Checking DOFS")
        self.Create_K()
        K = self.K
        hcount = 0
        for i in range(len(el)):
            if el[i].hybrid == 1:
                boun = np.hstack(el[i].BOUN)
                if self.ModelSet[0] == 'Planar':
                    for cdof in range(len(el[i].controlDOF)):            
                        if el[i].controlDOF[cdof] == 0:
                            
                            relevantDOFS = np.nonzero(self.Trt[self.hybridDOFS[cdof]])
                            relevantDOFS = relevantDOFS[0]
                            print(relevantDOFS)
                            controlneeded = 0
                            AnalyticalForces = 0
                            for m in range(len(relevantDOFS)):

                                if relevantDOFS[m]<=6*i+5 and relevantDOFS[m]>=6*i: #We are within hybrid element
                                    controlneeded+=1

                                elif el[int(np.floor(relevantDOFS[m]/6))].REL[int(relevantDOFS[m]%6)] == 0:
                                    print("???")
                                    pass
                                else:
                                    AnalyticalForces+=1

                            
                            if boun[cdof] == 1:
                                pass
                            elif (cdof == 0 or cdof == 3) and (boun[0] == 1 or boun[3] == 1):
                                print("Bound DOF")
                                pass
                            elif AnalyticalForces == 0:                            
                                print("Trivial DOF")
                                
                                self.K = np.delete(self.K,self.hybridDOFS[hcount*5+cdof],axis = 1)
                                self.K = np.delete(self.K,self.hybridDOFS[hcount*5+cdof],axis = 0)
                                self.Kh = np.delete(self.Kh,self.hybridDOFS[hcount*5+cdof],axis = 1)
                                self.Kh = np.delete(self.Kh,self.hybridDOFS[hcount*5+cdof],axis = 0)
                                self.M = np.delete(self.M,self.hybridDOFS[hcount*5+cdof],axis = 1)
                                self.M = np.delete(self.M,self.hybridDOFS[hcount*5+cdof],axis = 0)
                                self.P = np.delete(self.P,self.hybridDOFS[hcount*5+cdof])
                                self.FDOF = np.delete(self.FDOF,self.hybridDOFS[hcount*5+cdof])
                            elif controlneeded>0:
                                print(cdof)
                                print("ERROR: Not enough Controlled DOFS")
                            else:
                                print("Unknown Error, Model will not work as built")

                hcount+=1
                            
                            
            
        
        
            
    


        
        
    
    