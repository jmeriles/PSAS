#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:54:53 2020

@author: jim53
"""
import numpy as np
import sympy as sym
from scipy import integrate


class Element:
    
    node=[]
    nf =[]
    L=[]
    A=[]
    Ae = []
    BOUN=[]
    CON=[]
    ke=[]
    fe=[]
    REL=[]
    Iy=[]
    Iz=[]
    A=[]
    E=[]
    G=[]
    J=[]
    dx=[]
    dy=[]
    dz=[]
    id_v=[]
    v = []
    th = []
    Be = []
    hybrid = []
    controlDOF = []
    custom = 0
    
    def __init__(self,ElementType):
        self.ElementType=ElementType

    def get_Length(self):
    
        self.dx=self.node[1][0]-self.node[0][0]
        self.dy=self.node[1][1]-self.node[0][1]
        
        if self.ElementType[1]!='Planar':
            self.dz=self.node[1][2]-self.node[0][2]
        else:
            self.dz=0
            
        self.L=np.sqrt(self.dx*self.dx+self.dy*self.dy+self.dz*self.dz)
    
    def get_Area(self):
        if self.ElementType[0] =='Tri':
            self.Ae = (1/2)*((self.node[0][0])*(self.node[1][1]-self.node[2][1])+self.node[1][0]*(self.node[2][1]-self.node[0][1])\
                     +self.node[2][0]*(self.node[0][1]-self.node[1][1]))
        elif self.ElementType[0] == '4nQuad':
            x1 = self.node[0][0]
            x2 = self.node[1][0]
            x3 = self.node[2][0]
            x4 = self.node[3][0]
            y1 = self.node[0][1]
            y2 = self.node[1][1]
            y3 = self.node[2][1]
            y4 = self.node[3][1]
            self.Ae = (1/2)*((x1*y2+x2*y3+x3*y4+x4*y1)-(x2*y1+x3*y2+x4*y3+x1*y4))
                
    def Create_Rot(self):
        #     creates a rotation matrix that orients the element correctly and
        #     converts the forces from the local coordinate system to the global
        #     coordinate system
        
        #This shouldnt be necessary for the tri feas
        
        if self.ElementType[1]!='Planar':
            if self.dx==0 and self.dz==0:
                Rot=np.matrix([[0,1,0],[-1,0,0],[0,0,1]])
                Rot=np.transpose(Rot)
            else:
                lx=[self.dx,self.dy,self.dz]
                lx=lx/np.linalg.norm(lx)
                lz=np.cross(lx,[0,1,0])
                lz=lz/np.linalg.norm(lz)
                ly=np.cross(lx,lz) 
                ly=-ly/np.linalg.norm(ly)
                Rot=np.matrix([lx,ly,lz])
                Rot=np.transpose(Rot)
    #this is a rotation matrix I made manually, also works
        #         Rot=[dx,dy,dz;
        #        -dx*dy/sqrt(dx^2+dz^2),sqrt(dx^2+dz^2),-dz*dy/sqrt(dx^2+dz^2);
        #         -sqrt(dx^2+dy^2+dz^2)*dz/sqrt(dx^2+dz^2),0,sqrt(dx^2+dy^2+dz^2)*dx/sqrt(dx^2+dz^2)];
        #  
        #     Rot=Rot/L;
        else: 
            Rot=np.matrix([[self.dx/self.L,self.dy/self.L,0],\
                           [-self.dy/self.L,self.dx/self.L,0],\
                           [0,0,1]])
            Rot=np.transpose(Rot)
                
    #     Expands the Local to global rotation matrix to all twelve element
    #     forces
        r=0;
        self.ROT=Rot
        
        if self.ElementType[1]=='Planar':
            br=np.zeros((6,6))
        else:
            br=np.zeros((12,12))
        
        for n in range(int(len(br)/3)):
            for j in range(len(Rot)):
                for p in range(len(Rot)):
                    br[r+p,r+j]=Rot[p,j]
            r=r+j+1
            
        return [self.ROT,br]
    
        
    
    
    def create_Shape(self):
        if self.ElementType[0]=='Timoshenko':
            #The areas should be Ay and Az, but these depend on shape
            Phiy=12*self.E*self.Iz/(self.L**2*self.G*self.A)
            Phiz=12*self.E*self.Iy/(self.L**2*self.G*self.A)
        else:
            Phiy=0
            Phiz=0
        
        
        
        
        if (self.ElementType[0] == 'Timoshenko' or self.ElementType[0] == 'Bernoulli'):
            miuy=1/(1+Phiy)
            miuz=1/(1+Phiz)
            x=np.linspace(0,self.L,1000)
            Nu11=1-x/self.L
            Nu22=(1-miuy*(Phiy*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3))
            Nu62=(miuy*self.L/2)*((2+Phiy)*(x/self.L)-(4+Phiy)*(x/self.L)**2+2*(x/self.L)**3)
            Nu71=(x/self.L)
            Nu82=miuy*(Phiy*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3)
            Nu122=(miuy*self.L/2)*(-Phiy*(x/self.L)+(Phiy-2)*(x/self.L)**2+2*(x/self.L)**3)
            
        if self.ElementType[0] =='Tri':
            

            x11 = self.node[0][0]
            x12 = self.node[1][0]
            x13 = self.node[2][0]
            x21 = self.node[0][1]
            x22 = self.node[1][1]
            x23 = self.node[2][1]
            xcog = (x11+x12+x13)/3
            ycog = (x21+x22+x23)/3
            x11 = self.node[0][0]-xcog
            x12 = self.node[1][0]-xcog
            x13 = self.node[2][0]-xcog
            x21 = self.node[0][1]-ycog
            x22 = self.node[1][1]-ycog
            x23 = self.node[2][1]-ycog
            
            self.Be = 1/(2*self.A)*\
                 np.matrix([[x22-x23,0,x23-x21,0,x21-x22,0],\
                 [0,x13-x12,0,x11-x13,0,x12-x11],\
                 [x13-x12,x22-x23,x11-x13,x23-x21,x12-x11,x21-x22]])
                    
            N1 = x12*x23-x13*x22
                     
            N2 = x13*x21-x11*x23
                 
            N3 = x11*x22-x12*x21
            #none of this is necessary as per reader, but I wanted to practice integrating 


            
            self.Nu_int = 1/(2*self.Ae)*np.array([[N1,0,N2,0,N3,0],[0,N1,0,N2,0,N3]])
            
        if self.ElementType[0] =='4nQuad':
            x11 = self.node[0][0]
            x12 = self.node[1][0]
            x13 = self.node[2][0]
            x14 = self.node[3][0]
            x21 = self.node[0][1]
            x22 = self.node[1][1]
            x23 = self.node[2][1]
            x24 = self.node[3][1]
            
            xc = np.array([[x11,x21,x12,x22,x13,x23,x14,x24]])
            xc = xc.T
            
            e = sym.symbols('e')
            n = sym.symbols('n')
            
            N1 = (1/4)*(1-e)*(1-n)
            N2 = (1/4)*(1+e)*(1-n)
            N3 = (1/4)*(1+e)*(1+n)
            N4 = (1/4)*(1-e)*(1+n)
            
            dN1de = (1/4)*(-1)*(1-n)
            dN2de = (1/4)*(1)*(1-n)
            dN3de = (1/4)*(1)*(1+n)
            dN4de = (1/4)*(-1)*(1+n)
            dN1dn = (1/4)*(-1)*(1-e)
            dN2dn = (1/4)*(-1)*(1+e)
            dN3dn = (1/4)*(1)*(1+e)
            dN4dn = (1/4)*(1)*(1-e)
            
            self.Nu = np.array([[N1,0,N2,0,N3,0,N4,0],\
                          [0,N1,0,N2,0,N3,0,N4]])
                
            
            self.x = self.Nu @ xc
            
            dx1de = self.x[0][0].diff(e)
            dx1dn = self.x[0][0].diff(n)
            dx2de = self.x[1][0].diff(e)
            dx2dn = self.x[1][0].diff(n)
            
            self.J = np.array([[dx1de,dx2de],\
                          [dx1dn,dx2dn]])
            Jinv = abs(1/(dx1de*dx2dn-dx1dn*dx2de))*np.array([[dx2dn,-dx2de],
                                                              [-dx1dn,dx1de]])
            [dN1dx1,dN1dx2] = (Jinv)*np.array([[dN1de],[dN1dn]])
            [dN2dx1,dN2dx2] = (Jinv)*np.array([[dN2de],[dN2dn]]) 
            [dN3dx1,dN3dx2] = (Jinv)*np.array([[dN3de],[dN3dn]]) 
            [dN4dx1,dN4dx2] = (Jinv)*np.array([[dN4de],[dN4dn]]) 
            
            self.Be = np.array([[dN1dx1[0],0,dN2dx1[0],0,dN3dx1[0],0,dN4dx1[0],0],\
                                [0,dN1dx2[0],0,dN2dx2[0],0,dN3dx2[0],0,dN4dx2[0]],\
                                [dN1dx2[0],dN1dx1[0],dN2dx2[0],dN2dx1[0],dN3dx2[0],dN3dx1[0],dN4dx2[0],dN4dx1[0]]])
            
            
        if (self.ElementType[0] =='Timoshenko' or self.ElementType[0] =='Bernoulli'):
            if self.ElementType[1]!='Planar': 
                Nu33=(1-miuz*(Phiz*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3))
                Nu53=(-miuz*self.L/2)*((2+Phiz)*(x/self.L)-(4+Phiz)*(x/self.L)**2+2*(x/self.L)**3)
                Nu93=miuz*(Phiz*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3)
                Nu113=(-miuz*self.L/2)*(-Phiz*(x/self.L)+(Phiz-2)*(x/self.L)**2+2*(x/self.L)**3)
               
                
                self.Nu_int=np.matrix([[np.trapz(Nu11,x),0,0],\
                              [0,np.trapz(Nu22,x),0],\
                              [0,0,np.trapz(Nu33,x)],\
                              [0,0,0],\
                              [0,0,np.trapz(Nu53,x)],\
                              [0,np.trapz(Nu62,x),0],\
                              [np.trapz(Nu71,x),0,0],\
                              [0,np.trapz(Nu82,x),0],\
                              [0,0,np.trapz(Nu93,x)],\
                              [0,0,0],\
                              [0,0,np.trapz(Nu113,x)],\
                              [0,np.trapz(Nu122,x),0]])
                                  
                self.Nu_int=np.transpose(self.Nu_int)
                
            else:
                self.Nu_int=np.matrix([[np.trapz(Nu11,x),0],\
                              [0,np.trapz(Nu22,x)],\
                              [0,np.trapz(Nu62,x)],\
                              [np.trapz(Nu71,x),0],\
                              [0,np.trapz(Nu82,x)],\
                              [0,np.trapz(Nu122,x)]])
                                  
                self.Nu_int=np.transpose(self.Nu_int)
            
            Nt23=6*miuy*(-(x/self.L)+(x/self.L)**2)/self.L
            Nt63=1-miuy*((4+Phiy)*(x/self.L)-3*(x/self.L)**2)
            Nt83=6*miuy*((x/self.L)-(x/self.L)**2)/self.L
            Nt123=miuy*((Phiy-2)*(x/self.L)+3*(x/self.L)**2)
            
            if self.ElementType[1]!='Planar':
                Nt32=6*miuz*((x/self.L)-(x/self.L)**2)/self.L
                Nt41=1-(x/self.L)
                Nt52=1-miuz*((4+Phiz)*(x/self.L)-3*(x/self.L)**2)
                Nt92=6*miuz*(-(x/self.L)+(x/self.L)**2)/self.L
                Nt101=(x/self.L)
                Nt112=miuz*((Phiz-2)*(x/self.L)+3*(x/self.L)**2)
                
                
                self.Nt_int=np.matrix([[0,0,0],\
                              [0,0,np.trapz(Nt23,x)],\
                              [0,np.trapz(Nt32,x),0],\
                              [np.trapz(Nt41,x),0,0],
                              [0,np.trapz(Nt52,x),0],\
                              [0,0,np.trapz(Nt63,x)],\
                              [0,0,0],\
                              [0,0,np.trapz(Nt83,x)],\
                              [0,np.trapz(Nt92,x),0],\
                              [np.trapz(Nt101,x),0,0],\
                              [0,np.trapz(Nt112,x),0],\
                              [0,0,np.trapz(Nt123,x)]])
                self.Nt_int=np.transpose(self.Nt_int)
            else:
                self.Nt_int=np.matrix([[0],\
                              [np.trapz(Nt23,x)],\
                              [np.trapz(Nt63,x)],\
                              [0],\
                              [np.trapz(Nt83,x)],\
                              [np.trapz(Nt123,x)]])
                self.Nt_int=np.transpose(self.Nt_int)
                
               
            #Symbolic shape functions (For Later)
            x=sym.symbols('x')
            
            Nu11=1-x/self.L
            Nu22=(1-miuy*(Phiy*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3))
            Nu62=(miuy*self.L/2)*((2+Phiy)*(x/self.L)-(4+Phiy)*(x/self.L)**2+2*(x/self.L)**3)
            Nu71=(x/self.L)
            Nu82=miuy*(Phiy*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3)
            Nu122=(miuy*self.L/2)*(-Phiy*(x/self.L)+(Phiy-2)*(x/self.L)**2+2*(x/self.L)**3)
            
            if self.ElementType[1]!='Planar':
                Nu33=(1-miuz*(Phiz*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3))
                Nu53=(-miuz*self.L/2)*((2+Phiz)*(x/self.L)-(4+Phiz)*(x/self.L)**2+2*(x/self.L)**3)
                Nu93=miuz*(Phiz*(x/self.L)+3*(x/self.L)**2-2*(x/self.L)**3)
                Nu113=(-miuz*self.L/2)*(-Phiz*(x/self.L)+(Phiz-2)*(x/self.L)**2+2*(x/self.L)**3)
               
                
                self.Nu=np.matrix([[Nu11,0,0],\
                              [0,Nu22,0],\
                              [0,0,Nu33],\
                              [0,0,0],\
                              [0,0,Nu53],\
                              [0,Nu62,0],\
                              [Nu71,0,0],\
                              [0,Nu82,0],\
                              [0,0,Nu93],\
                              [0,0,0],\
                              [0,0,Nu113],\
                              [0,Nu122,0]])
                                  
                self.Nu=np.transpose(self.Nu)
                
            else:
                self.Nu=np.matrix([[Nu11,0],\
                              [0,Nu22],\
                              [0,Nu62],\
                              [Nu71,0],\
                              [0,Nu82],\
                              [0,Nu122]])
                                  
                self.Nu=np.transpose(self.Nu)
            
            Nt23=6*miuy*(-(x/self.L)+(x/self.L)**2)/self.L
            Nt63=1-miuy*((4+Phiy)*(x/self.L)-3*(x/self.L)**2)
            Nt83=6*miuy*((x/self.L)-(x/self.L)**2)/self.L
            Nt123=miuy*((Phiy-2)*(x/self.L)+3*(x/self.L)**2)
            
            if self.ElementType[1]!='Planar':
                Nt32=6*miuz*((x/self.L)-(x/self.L)**2)/self.L
                Nt41=1-(x/self.L)
                Nt52=1-miuz*((4+Phiz)*(x/self.L)-3*(x/self.L)**2)
                Nt92=6*miuz*(-(x/self.L)+(x/self.L)**2)/self.L
                Nt101=(x/self.L)
                Nt112=miuz*((Phiz-2)*(x/self.L)+3*(x/self.L)**2)
                
                
                self.Nt=np.matrix([[0,0,0],\
                              [0,0,Nt23],\
                              [0,Nt32,0],\
                              [Nt41,0,0],
                              [0,Nt52,0],\
                              [0,0,Nt63],\
                              [0,0,0],\
                              [0,0,Nt83],\
                              [0,Nt92,0],\
                              [Nt101,0,0],\
                              [0,Nt112,0],\
                              [0,0,Nt123]])
                self.Nt=np.transpose(self.Nt)
            else:
                self.Nt=np.matrix([[0],\
                              [Nt23],\
                              [Nt63],\
                              [0],\
                              [Nt83],\
                              [Nt123]])
                self.Nt=np.transpose(self.Nt)
        
    
    def create_LVec(self):
        x=sym.symbols('x')
        if self.ElementType[1]!='Planar':
            q=[self.w_x,self.w_y,self.w_z]
            m=[self.m_x,self.m_y,self.m_z]
        else:
            q=[self.w_x,self.w_y]
            m=[self.m_z]
        
        
        Nq=np.matmul(np.transpose(self.Nu_int),q)
        Ntt=np.matmul(np.transpose(self.Nt_int),m)
        
        self.p_w=Nq+Ntt
            
            
            
    def create_Q0(self):
        
        if self.ElementType[1]!='Planar':
            #Something here is wrong...
            self.eQ0=np.matrix([[self.q0,self.t0,-self.p_w[0,4],-self.p_w[0,10],-self.p_w[0,5],-self.p_w[0,11]]])
        else:
            self.eQ0=np.matrix([[self.q0,-self.p_w[0,2],-self.p_w[0,5]]])
                
    
    def release(self):
        L=self.L
        
        if self.ElementType[1]=='Planar':
            b=np.matrix([[-1,0,0],\
                          [0,1/L,1/L],\
                          [0,1,0],\
                          [1,0,0],\
                          [0,-1/L,-1/L],\
                          [0,0,1]])
            r=np.zeros(6)
            r=self.REL
        else: 
            b=np.matrix([[-1,0,0,0,0,0],\
             [0,0,0,0,1/L,1/L],\
             [0,0,-1/L,-1/L,0,0],\
             [0,-1,0,0,0,0],\
             [0,0,1,0,0,0],\
             [0,0,0,0,1,0],\
             [1,0,0,0,0,0],\
             [0,0,0,0,-1/L,-1/L],\
             [0,0,1/L,1/L,0,0],\
             [0,1,0,0,0,0],\
             [0,0,0,1,0,0],\
             [0,0,0,0,0,1]])
            r=np.zeros(12)
                    
                    
            r=self.REL
        
        zero_index=[]
        r_fixed=np.zeros(len(r))
        trigger=0
        for j in range(len(r)):
            if r[j]==0:
                zero_index.append(j)
                r_fixed[j]=np.inf
                
                trigger=1
            else:
                r_fixed[j]=r[j]
                if trigger==1:
                    trigger=1
                else:
                    trigger=0
        #Creates a vector for rearranging indices
        
        if self.ElementType[1]=='Planar':
            ind1=np.zeros(4)
            for k in range(len(zero_index)):
                if zero_index[k]==1:
                    ind1[0]=1
                if zero_index[k]==2:
                    ind1[1]=1
                if zero_index[k]==4:
                    ind1[2]=1
                if zero_index[k]==5:
                    ind1[3]=1
        else:
            ind2=np.zeros(4)
            ind3=np.zeros(4)            
            for k in range(len(zero_index)):
                if zero_index[k]==2:
                    ind2[0]=1
                if zero_index[k]==4:
                    ind2[1]=1
                if zero_index[k]==8:
                    ind2[2]=1
                if zero_index[k]==10:
                    ind2[3]=1
                if zero_index[k]==1:
                    ind3[0]=1
                if zero_index[k]==5:
                    ind3[1]=1
                if zero_index[k]==7:
                    ind3[2]=1
                if zero_index[k]==11:
                    ind3[3]=1
        
    
        
        if self.hybrid == 0:
            print(self.ke)
            self.ke_f=np.matmul(b,np.matmul(self.ke,np.transpose(b)))
        else:
            self.ke_f = self.ke

        
        #step 1, we set all the released indices to infinity for now
        
        
        
        #Reorder the stiffness matrix to three simpler matrices
        #r for restrained
        if self.ElementType[1]=='Planar':
            new_order=[0,3]
            Kr1 = [[self.ke_f[i,j] for j in new_order] for i in new_order]
            Kr1 = np.matrix(Kr1)
            new_order=[1,2,4,5]
            Kr3 = [[self.ke_f[i,j] for j in new_order] for i in new_order]
            Kr3 = np.matrix(Kr3)
            
            index=[0,3,1,2,4,5]
            self.T_i=np.zeros([6,6])
            #The T_i matrix is the transform matrix for any nonfixed springs
            #It only doesn't apply if we have a full release in which step 2 is 
            #necessary
            #Create fixed transform matrix
            
            for i in range (len(r_fixed)):
                for j in range(len(r_fixed)):
                    self.T_i[i,j]=self.ke_f[index[i],index[j]]/r_fixed[index[i]]
            
            self.T_i=self.T_i+np.identity(6)
            
            
            new_order=[0,3]
            T1_i = [[self.T_i[i,j] for j in new_order] for i in new_order]
            T1_i = np.matrix(T1_i)
            T1=np.linalg.inv(T1_i)
            new_order=[1,2,4,5]
            T3_i = [[self.T_i[i,j] for j in new_order] for i in new_order]
            T3_i = np.matrix(T3_i)
            T3=np.linalg.inv(T3_i)
            
            Tinf=np.bmat([[T1,np.zeros([2,4])],\
                      [np.zeros([4,2]),T3]])
            new_order=[0,2,3,1,4,5]
            Tinf = [[Tinf[i,j] for j in new_order] for i in new_order]
            
        else:
            new_order=[0,6,3,9]
            Kr1 = [[self.ke_f[i,j] for j in new_order] for i in new_order]
            Kr1 = np.matrix(Kr1)
            new_order=[2,4,8,10]
            Kr2 = [[self.ke_f[i,j] for j in new_order] for i in new_order]
            Kr2 = np.matrix(Kr2)
            new_order=[1,5,7,11]
            Kr3 = [[self.ke_f[i,j] for j in new_order] for i in new_order]
            Kr3 = np.matrix(Kr3)
            
            index=[0,6,3,9,2,4,8,10,1,5,7,11]
            self.T_i=np.zeros([12,12])
            #The T_i matrix is the transform matrix for any nonfixed springs
            #It only doesn't apply if we have a full release in which step 2 is 
            #necessary
            #Create fixed transform matrix
            
            for i in range (len(r_fixed)):
                for j in range(len(r_fixed)):
                    self.T_i[i,j]=self.ke_f[index[i],index[j]]/r_fixed[index[i]]
            
            self.T_i=self.T_i+np.identity(12)
            
            
            new_order=[0,1,2,3]
            T1_i = [[self.T_i[i,j] for j in new_order] for i in new_order]
            T1_i = np.matrix(T1_i)
            T1=np.linalg.inv(T1_i)
            new_order=[4,5,6,7]
            T2_i = [[self.T_i[i,j] for j in new_order] for i in new_order]
            T2_i = np.matrix(T2_i)
            T2=np.linalg.inv(T2_i)
            new_order=[8,9,10,11]
            T3_i = [[self.T_i[i,j] for j in new_order] for i in new_order]
            T3_i = np.matrix(T3_i)
            T3=np.linalg.inv(T3_i)
            
            Tinf=np.bmat([[T1,np.zeros([4,8])],\
                      [np.zeros([4,4]),T2,np.zeros([4,4])],\
                      [np.zeros([4,8]),T3]])
            new_order=[0,8,4,2,5,9,1,10,6,3,7,11]
            Tinf = [[Tinf[i,j] for j in new_order] for i in new_order]
        
        
        #Step 2, only needs to occur if we have full releases
        #f for free
        
        if self.ElementType[1]=='Planar':
            Kf3=np.matmul(Kr3,T3)
            
            if trigger==1:
                
                rdof3=[]
                cdof3=[]
                trigger3=0
                
                for m in range(4):
                    
                    if ind1[m]==1:
                        rdof3.append(m)
                        trigger3=1
                    else:
                        cdof3.append(m)
                        if trigger3==1:
                            trigger3=1
                        else:
                            trigger3=0
                            
                new_order3=rdof3+cdof3
            
                  
                if trigger3==1:
                    T_upright3=np.linalg.solve([[Kf3[i,j] for i in rdof3] for j in rdof3] ,  [[Kf3[i,j] for i in cdof3] for j in rdof3])
                    Tzero3=np.bmat([[np.matrix(np.zeros(((len(rdof3)),len(rdof3)))),-T_upright3],\
                                    [np.matrix(np.zeros((len(cdof3),len(rdof3)))),np.identity(len(cdof3))]])
                    
                    undo_order3=[new_order3.index(0),new_order3.index(1),new_order3.index(2),new_order3.index(3)]
                    Tzero3=[[Tzero3[j,i] for i in undo_order3] for j in undo_order3]
                else:
                    Tzero3=np.identity(4)
                    
                
                 
                
                Tzero1=np.identity(2)
                
                self.Tzero=np.bmat([[Tzero1,np.zeros([2,4])],\
                      [np.zeros([4,2]),Tzero3]])
                
                
            else:
                Tzero1=np.identity(2)
                Tzero3=np.identity(4)
                self.Tzero=np.identity(12)
            
            Kt1=np.matmul(Kr1,T1)
            Kt3=np.matmul(np.transpose(Tzero3),np.matmul(Kf3,Tzero3))
            
            Kt=np.bmat([[Kt1,np.zeros([2,4])],\
                        [np.zeros([4,2]),Kt3]])
                
            undo_order=[0,2,3,1,4,5]
            
            self.ke=np.matrix([[Kt[i,j] for i in undo_order] for j in undo_order])
            
            new_order=[0,3,1,2,4,5]
            self.p_w=np.matrix(self.p_w)
            self.p_w=np.matrix([self.p_w[0,i] for i in new_order])
            self.T=np.bmat([[np.matmul(T1,Tzero1),np.zeros([2,4])],\
                      [np.zeros([4,2]),np.matmul(T3,Tzero3)]])
                
            self.p_w=np.matmul(np.transpose(self.T),np.transpose(self.p_w))
            undo_order=[0,2,3,1,4,5]
            self.p_w=np.matrix([self.p_w[i,0] for i in undo_order])
            
            new_order=[0,3,1,2,4,5]
            self.me = [[self.me[i,j] for j in new_order] for i in new_order]
            self.me = np.matrix(self.me)
            self.me = np.matmul(np.transpose(self.T),np.matmul(self.me,self.T))
            
            undo_order=[0,2,3,1,4,5]
            self.me=np.matrix([[self.me[i,j] for i in undo_order] for j in undo_order])
        
           
            
            return [T1, T3,Tzero1, Tzero3]
        else:
            
            Kf2=np.matmul(Kr2,T2)
            Kf3=np.matmul(Kr3,T3)
                
            if trigger==1:
                
                
                rdof2=[]
                cdof2=[]
                rdof3=[]
                cdof3=[]
                #Rearrange the K matrices to put the released dofs in the upper left corner
                trigger2=0
                trigger3=0
                for m in range(4):
                    if ind2[m]==1:
                        rdof2.append(m)
                        trigger2=1
                    else:
                        cdof2.append(m)
                        if trigger2==1:
                            trigger2=1
                        else:
                            trigger2=0
                    if ind3[m]==1:
                        rdof3.append(m)
                        trigger3=1
                    else:
                        cdof3.append(m)
                        if trigger3==1:
                            trigger3=1
                        else:
                            trigger3=0
                            
                new_order2=rdof2+cdof2
                new_order3=rdof3+cdof3
                
                # if new_order2!=[0,1,2,3]:
                #     Kr2 = [[Kr2[i,j] for j in new_order2] for i in new_order2]
                #     Kr2 = np.matrix(Kr2)
                    
                # if new_order3!=[0,1,2,3]:
                #     Kr3 = [[Kr3[i,j] for j in new_order3] for i in new_order3]
                #     Kr3 = np.matrix(Kr3)
                
                #now use rearranged matrices to create the transformation matrix
                
                if trigger2==1:
                    T_upright2=np.linalg.solve([[Kf2[i,j] for i in rdof2] for j in rdof2] ,  [[Kf2[i,j] for i in cdof2] for j in rdof2])
                    Tzero2=np.bmat([[np.matrix(np.zeros(((len(rdof2)),len(rdof2)))),-T_upright2],\
                                    [np.matrix(np.zeros((len(cdof2),len(rdof2)))),np.identity(len(cdof2))]])
                    #rearrange back to original order
                    
                    undo_order2=[new_order2.index(0),new_order2.index(1),new_order2.index(2),new_order2.index(3)]
                    Tzero2=[[Tzero2[j,i] for i in undo_order2] for j in undo_order2]
                else:
                    Tzero2=np.identity(4)
                  
                if trigger3==1:
                    T_upright3=np.linalg.solve([[Kf3[i,j] for i in rdof3] for j in rdof3] ,  [[Kf3[i,j] for i in cdof3] for j in rdof3])
                    Tzero3=np.bmat([[np.matrix(np.zeros(((len(rdof3)),len(rdof3)))),-T_upright3],\
                                    [np.matrix(np.zeros((len(cdof3),len(rdof3)))),np.identity(len(cdof3))]])
                    
                    undo_order3=[new_order3.index(0),new_order3.index(1),new_order3.index(2),new_order3.index(3)]
                    Tzero3=[[Tzero3[j,i] for i in undo_order3] for j in undo_order3]
                else:
                    Tzero3=np.identity(4)
                    
                
                 
                
                Tzero1=np.identity(4)
                
                self.Tzero=np.bmat([[Tzero1,np.zeros([4,8])],\
                      [np.zeros([4,4]),Tzero2,np.zeros([4,4])],\
                      [np.zeros([4,8]),Tzero3]])
                # new_order=[11,7,5,1,10,8,4,2,9,3,6,0]
                    
                #     #should this be j,i?
                # Tzero = [[Tzero[j,i] for j in new_order] for i in new_order]
                # Tzero = np.matrix(Tzero)
                
            else:
                Tzero1=np.identity(4)
                Tzero2=np.identity(4)
                Tzero3=np.identity(4)
                self.Tzero=np.identity(12)
            
            Kt1=np.matmul(Kr1,T1)
            Kt2=np.matmul(np.transpose(Tzero2),np.matmul(Kf2,Tzero2))
            Kt3=np.matmul(np.transpose(Tzero3),np.matmul(Kf3,Tzero3))
            
            Kt=np.bmat([[Kt1,np.zeros([4,8])],\
                        [np.zeros([4,4]),Kt2,np.zeros([4,4])],\
                        [np.zeros([4,8]),Kt3]])
                
            undo_order=[0,8,4,2,5,9,1,10,6,3,7,11]
            
            self.ke=np.matrix([[Kt[i,j] for i in undo_order] for j in undo_order])
            
            new_order=[0,6,3,9,2,4,8,10,1,5,7,11]
            self.p_w=np.matrix(self.p_w)
            self.p_w=np.matrix([self.p_w[0,i] for i in new_order])
            self.T=np.bmat([[np.matmul(T1,Tzero1),np.zeros([4,8])],\
                      [np.zeros([4,4]),np.matmul(T2,Tzero2),np.zeros([4,4])],\
                      [np.zeros([4,8]),np.matmul(T3,Tzero3)]])
                
            self.p_w=np.matmul(np.transpose(self.T),np.transpose(self.p_w))
            undo_order=[0,8,4,2,5,9,1,10,6,3,7,11]
            self.p_w=np.matrix([self.p_w[i,0] for i in undo_order])
            
            new_order=[0,6,3,9,2,4,8,10,1,5,7,11]
            self.me = [[self.me[i,j] for j in new_order] for i in new_order]
            self.me = np.matrix(self.me)
            self.me = np.matmul(np.transpose(self.T),np.matmul(self.me,self.T))
            
            undo_order=[0,8,4,2,5,9,1,10,6,3,7,11]
            self.me=np.matrix([[self.me[i,j] for i in undo_order] for j in undo_order])
        
           
            
            return [T1, T2, T3,Tzero1,Tzero2,Tzero3]
        