#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 15:34:13 2020

@author: jim53
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D 

class PlotModel:
    
    def __init__(self):
        pass

    def Struct(self,Model,el):
        #Plotting
        #Plot Model
        #Take element end nodes and plot them in 3d space. Include BC and Loading
        #Vectors
        if Model.ModelSet[0]=='Planar' and Model.ModelSet[1] == 'line':
            x=[]
            y=[]
            xh = []
            yh = []
            
            maxL=0
            for i in range(Model.numel):
                if el[i].hybrid == 0:
                    x.append(el[i].node[0][0])
                    x.append(el[i].node[1][0])
                    y.append(el[i].node[0][1])
                    y.append(el[i].node[1][1])
                else:
                    xh.append(el[i].node[0][0])
                    xh.append(el[i].node[1][0])
                    yh.append(el[i].node[0][1])
                    yh.append(el[i].node[1][1])

                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
                
            
            #find maxL element length for scaling
        
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            if Model.hybrid == 1:
                ax.plot(xh,yh,'g--')
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
        
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            # ax.set_xlim3d(-5,5)
            # ax.set_ylim3d(-5,5)
        

            #Find how many BCs we have and try to represent them in the plot
            #straight lines generated for the translational fixities and circles 
            #for rotational fixities
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')
                        
                        if j==2:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r',linewidth=.75)
                    
            #Plot Release
            for i in range(Model.numel):
                for j in range(len(el[i].REL)):
                    if el[i].REL[j] != np.inf:
                        if j == 0:
                            x=np.matrix([el[i].L/2-el[i].L/10,el[i].L/2-el[i].L/10])
                            y=np.matrix([el[i].L/50,-el[i].L/50])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                            coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()
                            
                            if el[i].REL[j] == 0:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g')
                            else:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g',linestyle="--")
                            ax.add_line(line1)
                        if j == 3:
                            x=np.matrix([el[i].L/2+el[i].L/10,el[i].L/2+el[i].L/10])
                            y=np.matrix([el[i].L/50,-el[i].L/50])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                            coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()
                            
                            if el[i].REL[j] == 0:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g')
                            else:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g',linestyle="--")

                            ax.add_line(line1)
                            
                        if j == 1:
                            x=np.matrix([el[i].L/2-el[i].L/50,el[i].L/2+el[i].L/50])
                            y=np.matrix([el[i].L/50,el[i].L/50])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                            coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()
                            
                            if el[i].REL[j] == 0:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g')
                            else:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g',linestyle="--")
                            ax.add_line(line1)
                        if j == 4:
                            x=np.matrix([el[i].L/2-el[i].L/50,el[i].L/2+el[i].L/50])
                            y=np.matrix([-el[i].L/50,-el[i].L/50])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                            coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()
                            
                            if el[i].REL[j] == 0:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g')
                            else:
                                line1 = Line2D(coord_g[0], coord_g[1],color = 'g',linestyle="--")

                            ax.add_line(line1)
                            
                        if j == 2:
                            x=np.matrix([el[i].L/20])
                            y=np.matrix([0])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]

                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()

                            if el[i].REL[j] == 0:
                                circle1 = plt.Circle(tuple(np.hstack(coord_g)), .1, color='g',fill = 1)
                            else:
                                circle1 = plt.Circle(tuple(np.hstack(coord_g)), .1, color='g',fill = 0)
                            ax.add_patch(circle1)
                        if j == 5:
                            x=np.matrix([el[i].L-el[i].L/20])
                            y=np.matrix([0])
                    
                            coord_l=np.concatenate((x,y))
                            coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                            
                            node_i=np.array(el[i].node[0])
                            
                            coord_g[:,0]=coord_g[:,0]+np.c_[node_i]

                            
                            coord_l=coord_l.tolist()
                            coord_g=coord_g.tolist()
                            
                            if el[i].REL[j] == 0:
                                circle1 = plt.Circle(tuple(np.hstack(coord_g)), .1, color='g',fill = 1)
                            else:
                                circle1 = plt.Circle(tuple(np.hstack(coord_g)), .1, color='g',fill = 0)
                            ax.add_patch(circle1)
                            
            
            #Plots nodal loads 
            for i in range(len(Model.NODE)):
                for j in range(len(Model.LOAD[0])):
                    if Model.LOAD[i][j]!=0:
                        if j==0:
                            if Model.LOAD[i][j]>0:
                                Lx=np.linspace(0,maxL/10,2)+Model.NODE[i][0]
                            if Model.LOAD[i][j]<0:
                                Lx=np.linspace(0,maxL/10,2)+Model.NODE[i][0]
                            Ly=np.zeros(2)+Model.NODE[i][1]
                            
                            plt.arrow(Lx[0],Ly[0],Lx[1]-Lx[0],Ly[1]-Ly[0],width=.025)
                            ax.text(Lx[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                        
                        if j==1:
                            Lx=np.zeros(2)+Model.NODE[i][0]
                            if Model.LOAD[i][j]>0:
                                Ly=np.linspace(0,maxL/10,2)+Model.NODE[i][1]
                            if Model.LOAD[i][j]<0:
                                Ly=np.linspace(0,-maxL/10,2)+Model.NODE[i][1]
                            
                            plt.arrow(Lx[0],Ly[0],Lx[1]-Lx[0],Ly[1]-Ly[0],width=.025)
                            ax.text(Lx[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                       
                        if j==2:
                            if Model.LOAD[i][j]>0:
                                theta=np.linspace(0,1.9*np.pi,72)
                            if Model.LOAD[i][j]<0:
                                theta=np.linspace(0,-1.9*np.pi,72)
                            
                            Lx=maxL/20*np.cos(theta)+Model.NODE[i][0]
                            Ly=maxL/20*np.sin(theta)+Model.NODE[i][1]
                            
                            ax.plot(Lx,Ly, color='b',linewidth=.75)
                            
                            L_ax=[Lx[-2],Lx[-1]]
                            L_ay=[Ly[-2],Ly[-1]]
                            
                            plt.arrow(L_ax[0],L_ay[0],L_ax[1]-L_ax[0],L_ay[1]-L_ay[0],width=.025)
                            ax.text(Lx[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                            
                            
                            
              #Plot Distributed load
             
            #perp=[]
            for i in range(Model.numel):
                
                if el[i].w_y!=0:
                    x=np.matrix([0,el[i].L/4,el[i].L/2,3*el[i].L/4,el[i].L])
                    y=np.matrix([maxL/20,maxL/20,maxL/20,maxL/20,maxL/20])
            
                    coord_l=np.concatenate((x,y))
                    coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                    
                    node_i=np.array(el[i].node[0])
                    
                    coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                    coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                    coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                    coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                    coord_g[:,4]=coord_g[:,4]+np.c_[node_i]
                    
                    coord_l=coord_l.tolist()
                    coord_g=coord_g.tolist()
                    ax.plot(coord_g[0],coord_g[1], color='b',linewidth=.75)
                    ax.text(coord_g[0][2],coord_g[1][2],str('{:.3f}'.format(el[i].w_y)))
                    
                    for j in range(len(coord_l[1])):
                        
                        if el[i].w_y<0:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([y[0,j],0])
                        else:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([0,y[0,j]])
                        
                        arrow_coord_l=np.concatenate((arrow_x,arrow_y))
                        arrow_coord_g=np.matmul(el[i].ROT[0:2,0:2],arrow_coord_l)
                        
                        arrow_coord_g[:,0]=arrow_coord_g[:,0]+np.c_[node_i]
                        arrow_coord_g[:,1]=arrow_coord_g[:,1]+np.c_[node_i]
                        
                        arrow_coord_g=arrow_coord_g.tolist()
                        ax.plot(arrow_coord_g[0],arrow_coord_g[1],color='b',linewidth=.75)
                        plt.arrow(arrow_coord_g[0][0],arrow_coord_g[1][0],arrow_coord_g[0][1]-arrow_coord_g[0][0],arrow_coord_g[1][1]-arrow_coord_g[1][0],width=.025)
                        
            plt.show() 
        elif Model.ModelSet[0] == '3D' and Model.ModelSet[1] == 'line':
        
            x=[]
            y=[]
            z=[]
            xh=[]
            yh=[]
            zh=[]
            maxL=0
            for i in range(Model.numel):
                if el[i].hybrid ==0:
                    x.append(el[i].node[0][0])
                    x.append(el[i].node[1][0])
                    y.append(el[i].node[0][1])
                    y.append(el[i].node[1][1])
                    z.append(el[i].node[0][2])
                    z.append(el[i].node[1][2])
                else:
                    xh.append(el[i].node[0][0])
                    xh.append(el[i].node[1][0])
                    yh.append(el[i].node[0][1])
                    yh.append(el[i].node[1][1])
                    zh.append(el[i].node[0][2])
                    zh.append(el[i].node[1][2])
                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
                
            
            #find maxL element length for scaling
        
            
            from matplotlib.patches import FancyArrowPatch
            from mpl_toolkits.mplot3d import proj3d
            
            
            #Creates the arrows for plotting forces
            class Arrow3D(FancyArrowPatch):
                def __init__(self, xs, ys, zs, *args, **kwargs):
                    FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
                    self._verts3d = xs, ys, zs
            
                def draw(self, renderer):
                    xs3d, ys3d, zs3d = self._verts3d
                    xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
                    self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
                    FancyArrowPatch.draw(self, renderer)
                    
                    
    
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            if Model.hybrid == 1:
                ax.plot(xh,zh,yh,'g--')
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            # ax.set_xlim3d(-5,5)
            # ax.set_ylim3d(-5,5)
            
            
            #Find how many BCs we have and try to represent them in the plot
            #straight lines generated for the translational fixities and circles 
            #for rotational fixities
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            BCz=np.zeros(2)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            BCz=np.zeros(10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==2:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.zeros(10)+Model.NODE[i][1]
                            BCz=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=np.zeros(72)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=maxL/40*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==4:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=np.zeros(72)+Model.NODE[i][1]
                            BCz=maxL/40*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                    
            
            #Plots nodal loads but needs a better way to annotate the loads
            for i in range(len(Model.NODE)):
                for j in range(len(Model.LOAD[0])):
                    if Model.LOAD[i][j]!=0:
                        if j==0:
                            if Model.LOAD[i][j]>0:
                                Lx=np.linspace(0,maxL/10,2)+Model.NODE[i][0]
                            if Model.LOAD[i][j]<0:
                                Lx=np.linspace(0,maxL/10,2)+Model.NODE[i][0]
                            Ly=np.zeros(2)+Model.NODE[i][1]
                            Lz=np.zeros(2)+Model.NODE[i][2]
                            
                            a = Arrow3D(Lx,Lz,Ly,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],\
                                    Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                        
                        if j==1:
                            Lx=np.zeros(2)+Model.NODE[i][0]
                            if Model.LOAD[i][j]>0:
                                Ly=np.linspace(0,maxL/10,2)+Model.NODE[i][1]
                            if Model.LOAD[i][j]<0:
                                Ly=np.linspace(0,-maxL/10,2)+Model.NODE[i][1]
                            Lz=np.zeros(2)+Model.NODE[i][2]
                            
                            a = Arrow3D(Lx,Lz,Ly,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                        if j==2:
                            Lx=np.zeros(2)+Model.NODE[i][0]
                            Ly=np.zeros(2)+Model.NODE[i][1]
                            if Model.LOAD[i][j]>0:
                                Lz=np.linspace(0,maxL/10,2)+Model.NODE[i][2]
                            if Model.LOAD[i][j]<0:
                                Lz=np.linspace(0,-maxL/10,2)+Model.NODE[i][2]
                            
                            a = Arrow3D(Lx,Lz,Ly,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                        if j==3:
                            if Model.LOAD[i][j]>0:
                                theta=np.linspace(0,1.9*np.pi,72)
                            if Model.LOAD[i][j]<0:
                                theta=np.linspace(0,-1.9*np.pi,72)
                            Lx=np.zeros(72)+Model.NODE[i][0]
                            Ly=maxL/20*np.sin(theta)+Model.NODE[i][1]
                            Lz=maxL/20*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(Lx,Lz,Ly, color='b',linewidth=.75)
                            
                            L_ax=[Lx[-1],Lx[-1]]
                            L_ay=[Ly[-2],Ly[-1]]
                            L_az=[Lz[-2],Lz[-1]]
                            a = Arrow3D(L_ax,L_az,L_ay,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
            
                        if j==4:
                            if Model.LOAD[i][j]>0:
                                theta=np.linspace(0,1.9*np.pi,72)
                            if Model.LOAD[i][j]<0:
                                theta=np.linspace(0,-1.9*np.pi,72)
                            
                            Lx=maxL/20*np.cos(theta)+Model.NODE[i][0]
                            Ly=np.zeros(72)+Model.NODE[i][1]
                            Lz=maxL/20*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(Lx,Lz,Ly, color='b',linewidth=.75)
                            
                            L_ax=[Lx[-2],Lx[-1]]
                            L_ay=[Ly[-1],Ly[-1]]
                            L_az=[Lz[-2],Lz[-1]]
                            a = Arrow3D(L_ax,L_az,L_ay,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                            
                            
                        if j==5:
                            if Model.LOAD[i][j]>0:
                                theta=np.linspace(0,1.9*np.pi,72)
                            if Model.LOAD[i][j]<0:
                                theta=np.linspace(0,-1.9*np.pi,72)
                            
                            Lx=maxL/20*np.cos(theta)+Model.NODE[i][0]
                            Ly=maxL/20*np.sin(theta)+Model.NODE[i][1]
                            Lz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(Lx,Lz,Ly, color='b',linewidth=.75)
                            
                            L_ax=[Lx[-2],Lx[-1]]
                            L_ay=[Ly[-2],Ly[-1]]
                            L_az=[Lz[-1],Lz[-1]]
                            a = Arrow3D(L_ax,L_az,L_ay,mutation_scale=5,\
                                        lw=1, arrowstyle="-|>", color="b")
                            ax.add_artist(a)
                            ax.text(Lx[-1],Lz[-1],Ly[-1],str('{:.3f}'.format(Model.LOAD[i][j])))
                            
                            
             #Plot Distributed load
             
            #perp=[]
            for i in range(Model.numel):
                if el[i].w_z!=0:
                    x=np.matrix([0,el[i].L/4,el[i].L/2,3*el[i].L/4,el[i].L])
                    y=np.matrix([0,0,0,0,0])
                    z=np.matrix([maxL/20,maxL/20,maxL/20,maxL/20,maxL/20])
            
                    coord_l=np.concatenate((x,y,z))
                    coord_g=np.matmul(el[i].ROT,coord_l)
                    
                    node_i=np.array(el[i].node[0])
                    
                    coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                    coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                    coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                    coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                    coord_g[:,4]=coord_g[:,4]+np.c_[node_i]
                    
                    coord_l=coord_l.tolist()
                    coord_g=coord_g.tolist()
                    ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                    ax.text(coord_g[0][2],coord_g[2][2]\
                            ,coord_g[1][2],str('{:.3f}'.format(el[i].w_z)))
                    
                    for j in range(len(coord_l[1])):
                        
                        if el[i].w_z<0:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([y[0,j],0])
                            arrow_z=np.matrix([z[0,j],0])
                        else:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([0,y[0,j]])
                            arrow_z=np.matrix([0,z[0,j]])
                        
                        arrow_coord_l=np.concatenate((arrow_x,arrow_y,arrow_z))
                        arrow_coord_g=np.matmul(el[i].ROT,arrow_coord_l)
                        
                        arrow_coord_g[:,0]=arrow_coord_g[:,0]+np.c_[node_i]
                        arrow_coord_g[:,1]=arrow_coord_g[:,1]+np.c_[node_i]
                        
                        
                        arrow_coord_g=arrow_coord_g.tolist()
                        ax.plot(arrow_coord_g[0],arrow_coord_g[2],arrow_coord_g[1],\
                                color='b',linewidth=.75)
                        a = Arrow3D(arrow_coord_g[0],arrow_coord_g[2],arrow_coord_g[1],\
                                    mutation_scale=5,lw=1, arrowstyle="-|>", color="b")
                        ax.add_artist(a)
                if el[i].w_y!=0:
                    x=np.matrix([0,el[i].L/4,el[i].L/2,3*el[i].L/4,el[i].L])
                    y=np.matrix([maxL/20,maxL/20,maxL/20,maxL/20,maxL/20])
                    z=np.matrix([0,0,0,0,0])
            
                    coord_l=np.concatenate((x,y,z))
                    coord_g=np.matmul(el[i].ROT,coord_l)
                    
                    node_i=np.array(el[i].node[0])
                    
                    coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                    coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                    coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                    coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                    coord_g[:,4]=coord_g[:,4]+np.c_[node_i]
                    
                    coord_l=coord_l.tolist()
                    coord_g=coord_g.tolist()
                    ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                    ax.text(coord_g[0][2],coord_g[2][2]\
                            ,coord_g[1][2],str('{:.3f}'.format(el[i].w_y)))
                    
                    for j in range(len(coord_l[1])):
                        
                        if el[i].w_y<0:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([y[0,j],0])
                            arrow_z=np.matrix([z[0,j],0])
                        else:
                            arrow_x=np.matrix([x[0,j],x[0,j]])
                            arrow_y=np.matrix([0,y[0,j]])
                            arrow_z=np.matrix([0,z[0,j]])
                        
                        arrow_coord_l=np.concatenate((arrow_x,arrow_y,arrow_z))
                        arrow_coord_g=np.matmul(el[i].ROT,arrow_coord_l)
                        
                        arrow_coord_g[:,0]=arrow_coord_g[:,0]+np.c_[node_i]
                        arrow_coord_g[:,1]=arrow_coord_g[:,1]+np.c_[node_i]
                        
                        arrow_coord_g=arrow_coord_g.tolist()
                        ax.plot(arrow_coord_g[0],arrow_coord_g[2],arrow_coord_g[1],\
                                color='b',linewidth=.75)
                        a = Arrow3D(arrow_coord_g[0],arrow_coord_g[2],arrow_coord_g[1],\
                                    mutation_scale=5,lw=1, arrowstyle="-|>", color="b")
                        ax.add_artist(a)
                        
                    
                        
                        
                        #This checks if arrows are perpendicular... they are but they dont
                        #look like it
                        # a_vec=[arrow_coord_g[0][1]-arrow_coord_g[0][0],\
                        #        arrow_coord_g[1][1]-arrow_coord_g[1][0],\
                        #        arrow_coord_g[2][1]-arrow_coord_g[2][0]]
                        # e_vec=[coord_g[0][1]-coord_g[0][0],\
                        #        coord_g[1][1]-coord_g[1][0],\
                        #        coord_g[2][1]-coord_g[2][0]]
                        # if np.dot(a_vec,e_vec)<.00000001:
                        #     perp.append(1)
                        # else:
                        #     perp.append(0)
            plt.show()       
            ax.set_ylim3d(-5,0)
        elif Model.ModelSet[1] == 'FEA_Tri' :
            x=[]
            y=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                x.append(el[i].node[2][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                y.append(el[i].node[2][1])

                #find maxL element length for scaling
                # if el[i].L>maxL:
                #     maxL=el[i].L
                
            
            #find maxL element length for scaling
        
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
        
        elif Model.ModelSet[1] == 'FEA_4nQuad':
            fig=plt.figure()
            ax = fig.add_subplot(111)
            
            for i in range(Model.numel):
                x = []
                y = []
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                x.append(el[i].node[2][0])
                x.append(el[i].node[3][0])
                x.append(el[i].node[0][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                y.append(el[i].node[2][1])
                y.append(el[i].node[3][1])
                y.append(el[i].node[0][1])
                ax.plot(x,y,'b')
            
        
    def DefShape(self,Model,el):
        
        #Plot Deformed Shape
        
        if Model.ModelSet[0]=='Planar' and Model.ModelSet[1] =='line':
            #Extrapolate and plot Deformed Shape
            #Start with axial deformation
            
            x=[]
            y=[]
            
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])

            
            
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            
            for i in range(Model.numel):
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                Uel=np.zeros((1,Model.nf))
                
                for j in range(len(el_boun)):
                    if el_boun[j]==0:
                        Uel[0,j]=Model.U[Model.FDOF.index(el[i].id_v[j]-1)]
                    else:
                        Uel[0,j]=0
                
                #This should be a function (Create br)
                ibr=np.zeros((6,6))
                r=0
                for n in range(2):
                    for j in range(len(el[i].ROT)):
                        for p in range(len(el[i].ROT)):
                            ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                    r=r+j+1
                
                el_U_l=np.matmul(ibr,np.transpose(Uel))
                
                #Start with axial
                dL=el_U_l[3]-el_U_l[0]
                
                
                #Create Shape Functions
                xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[3],20))
                #x=np.transpose(np.linspace(0,el[i].L,20))
                x=x.squeeze()
                xl=np.matrix(xl)
                x=np.matrix(x)
                l_y=[]
                
                for n in range(x.size):
                    N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                    N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                    N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                    N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                    
                    N=np.array([N1,N2,N3,N4])
                    Uy=np.array([el_U_l[1],el_U_l[2],el_U_l[4],el_U_l[5]])
                    
                    l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                    
                    
                l_y=np.matrix(l_y)
                
                def_l=np.concatenate((x,l_y))
                
                def_g=np.matmul(el[i].ROT[0:2,0:2],def_l)
                def_g=def_g.tolist()
                
                def_g[0]=[x+el[i].node[0][0] for x in def_g[0]]
                def_g[1]=[x+el[i].node[0][1] for x in def_g[1]]
                
                ax.plot(def_g[0],def_g[1],color='r',linewidth=.75)
                
        elif Model.ModelSet[0]=='3D' and Model.ModelSet[1] =='line':
            #Extrapolate and plot Deformed Shape
            #Start with axial deformation
            
            x=[]
            y=[]
            z=[]
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                z.append(el[i].node[0][2])
                z.append(el[i].node[1][2])
            
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            for i in range(Model.numel):
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                Uel=np.zeros((1,Model.nf))
                for j in range(len(el_boun)):
                    if el_boun[j]==0:
                        Uel[0,j]=Model.U[Model.FDOF.index(el[i].id_v[j]-1)]
                    else:
                        Uel[0,j]=0
                
                #This should be a function (Create br)
                ibr=np.zeros((12,12))
                r=0
                for n in range(4):
                    for j in range(len(el[i].ROT)):
                        for p in range(len(el[i].ROT)):
                            ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                    r=r+j+1
                
                el_U_l=np.matmul(ibr,np.transpose(Uel))
                
                #Start with axial
                dL=el_U_l[6]-el_U_l[0]
                
                #Now work on z axis
                #Create Shape Functions
                xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[6],20))
                x=x.squeeze()
                x=np.matrix(x)
                xl=np.matrix(xl)
                l_z=[]
                l_y=[]
                
                for n in range(x.size):
                    N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                    N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                    N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                    N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                    
                    N=np.array([N1,N2,N3,N4])
                    Uz=np.array([el_U_l[2],el_U_l[4],el_U_l[8],el_U_l[10]])
                    Uy=np.array([el_U_l[1],el_U_l[5],el_U_l[7],el_U_l[11]])
                    
                    l_z.append(np.matmul(np.transpose(N),Uz)[0,0])
                    l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                    
                    
                l_z=np.matrix(l_z)
                l_y=np.matrix(l_y)
                
                def_l=np.concatenate((x,l_y,l_z))
                
                def_g=np.matmul(el[i].ROT,def_l)
                def_g=def_g.tolist()
                
                def_g[0]=[x+el[i].node[0][0] for x in def_g[0]]
                def_g[1]=[x+el[i].node[0][1] for x in def_g[1]]
                def_g[2]=[x+el[i].node[0][2] for x in def_g[2]]
                
                ax.plot(def_g[0],def_g[2],def_g[1],color='r',linewidth=.75)
        elif Model.ModelSet[1] == 'FEA_Tri':
            
            #Extrapolate and plot Deformed Shape
            #Start with axial deformation

            x=[]
            y=[]
            xo = []
            yo = []
            defx = []
            defy = []
            
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                x.append(el[i].node[2][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                y.append(el[i].node[2][1])

            
            
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            
            for i in range(Model.numel):
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]+el[i].BOUN[2]
                Uel=np.zeros((1,Model.nf))
                
                count = 0;
                for j in range(len(el_boun)):
    
                    if el_boun[j]==0 and count==0:
                        defx.append(float(Model.U[Model.FDOF.index(el[i].id_v[j]-1)]))
                        count = count+1
                    elif el_boun[j]==0 and count == 1:
                        defy.append(float(Model.U[Model.FDOF.index(el[i].id_v[j]-1)]))
                        count = 0
                    elif el_boun[j]==1 and count == 0:
                       defx.append(0)
                       count = count+1
                    elif el_boun[j]==1 and count == 1:  
                        defy.append(0)
                        count = 0
            defx = np.array(x)+np.array(defx)
            defy = np.array(y)+np.array(defy)
            
            print(defx)
            print(defy)
            ax.plot(defx.tolist(),defy.tolist(),'r')
        elif Model.ModelSet[1] == 'FEA_4nQuad':
            

            fig=plt.figure()
            ax = fig.add_subplot(111)
            
            for i in range(Model.numel):
                x = []
                y = []
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                x.append(el[i].node[2][0])
                x.append(el[i].node[3][0])
                x.append(el[i].node[0][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                y.append(el[i].node[2][1])
                y.append(el[i].node[3][1])
                y.append(el[i].node[0][1])
                ax.plot(x,y,'b')

            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            
            for i in range(Model.numel):
                xo = []
                yo = []
                defx = []
                defy = []
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]+el[i].BOUN[2]+el[i].BOUN[3]
                Uel=np.zeros((1,Model.nf))
                
                xo.append(el[i].node[0][0])
                xo.append(el[i].node[1][0])
                xo.append(el[i].node[2][0])
                xo.append(el[i].node[3][0])
                xo.append(el[i].node[0][0])
                yo.append(el[i].node[0][1])
                yo.append(el[i].node[1][1])
                yo.append(el[i].node[2][1])
                yo.append(el[i].node[3][1])
                yo.append(el[i].node[0][1])
                count = 0;
                for j in range(len(el_boun)):
    
                    if el_boun[j]==0 and count==0:
                        defx.append(float(Model.U[Model.FDOF.index(el[i].id_v[j]-1)]))
                        count = count+1
                    elif el_boun[j]==0 and count == 1:
                        defy.append(float(Model.U[Model.FDOF.index(el[i].id_v[j]-1)]))
                        count = 0
                    elif el_boun[j]==1 and count == 0:
                       defx.append(0)
                       count = count+1
                    elif el_boun[j]==1 and count == 1:  
                        defy.append(0)
                        count = 0
                defx.append(defx[0])
                defy.append(defy[0])
                defx = np.array(xo)+np.array(defx)
                defy = np.array(yo)+np.array(defy)
                ax.plot(defx,defy,'r')

            

                    
                
                

               
    def DynDefShape(self,Model,el,k):
        
        #this function takes a k value, and so should be able to get the deformed
        #shape at a specific time, but is not done yet
        if Model.ModelSet[0]=='Planar':
            #Plot deformed shape for a specific time
            #Extrapolate and plot Deformed Shape
            
            
            fig=plt.figure()
            ax = fig.add_subplot(111)
          
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            
            for i in range(Model.numel):
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                Uel=np.zeros((1,Model.nf))
                for j in range(len(el_boun)):
                    if el_boun[j]==0:
                        Uel[0,j]=Model.Udy[Model.FDOF.index(el[i].id_v[j]-1),k]
                    else:
                        Uel[0,j]=0
                
                #This should be a function (Create br)
                ibr=np.zeros((6,6))
                r=0
                for n in range(2):
                    for j in range(len(el[i].ROT)):
                        for p in range(len(el[i].ROT)):
                            ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                    r=r+j+1
                
                el_U_l=np.matmul(ibr,np.transpose(Uel))
                
                #Start with axial
                dL=el_U_l[3]-el_U_l[0]
                
                #Now work on z axis
                
                #Create Shape Functions
                xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[3],20))
                x=x.squeeze()
                xl=np.matrix(xl)
                x=np.matrix(x)
                l_y=[]
                
                for n in range(x.size):
                    N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                    N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                    N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                    N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                    
                    N=np.array([N1,N2,N3,N4])
                    Uy=np.array([el_U_l[1],el_U_l[2],el_U_l[4],el_U_l[5]])
                    
                    l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                    
                    
                l_y=np.matrix(l_y)
                
                def_l=np.concatenate((x,l_y))
                
                def_g=np.matmul(el[i].ROT[0:2,0:2],def_l)
                def_g=def_g.tolist()
                
                def_g[0]=[x+el[i].node[0][0] for x in def_g[0]]
                def_g[1]=[x+el[i].node[0][1] for x in def_g[1]]
                
                ax.plot(def_g[0],def_g[1],color='r',linewidth=.75)
        else:
            #Plot deformed shape for a specific time
            #Extrapolate and plot Deformed Shape
            
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
          
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            for i in range(Model.numel):
                el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                Uel=np.zeros((1,Model.nf))
                for j in range(len(el_boun)):
                    if el_boun[j]==0:
                        Uel[0,j]=Model.Udy[Model.FDOF.index(el[i].id_v[j]-1),k]
                    else:
                        Uel[0,j]=0
                
                #This should be a function (Create br)
                ibr=np.zeros((12,12))
                r=0
                for n in range(4):
                    for j in range(len(el[i].ROT)):
                        for p in range(len(el[i].ROT)):
                            ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                    r=r+j+1
                
                el_U_l=np.matmul(ibr,np.transpose(Uel))
                
                #Start with axial
                dL=el_U_l[6]-el_U_l[0]
                
                #Now work on z axis
                
                 #Create Shape Functions
                xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[6],20))
                x=x.squeeze()
                x=np.matrix(x)
                xl=np.matrix(xl)
                l_z=[]
                l_y=[]
                
                for n in range(x.size):
                    N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                    N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                    N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                    N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                    
                    N=np.array([N1,N2,N3,N4])
                    Uz=np.array([el_U_l[2],el_U_l[4],el_U_l[8],el_U_l[10]])
                    Uy=np.array([el_U_l[1],el_U_l[5],el_U_l[7],el_U_l[11]])
                    
                    l_z.append(np.matmul(np.transpose(N),Uz)[0,0])
                    l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                    
                    
                l_z=np.matrix(l_z)
                l_y=np.matrix(l_y)
                
                def_l=np.concatenate((x,l_y,l_z))
                
                def_g=np.matmul(el[i].ROT,def_l)
                def_g=def_g.tolist()
                
                def_g[0]=[x+el[i].node[0][0] for x in def_g[0]]
                def_g[1]=[x+el[i].node[0][1] for x in def_g[1]]
                def_g[2]=[x+el[i].node[0][2] for x in def_g[2]]
                
                ax.plot(def_g[0],def_g[2],def_g[1],color='r',linewidth=.75)
            

    def AnimateDefShape(self,Model,el,scale):
        from matplotlib.animation import FuncAnimation
        
        if Model.ModelSet[0]=='Planar':
            fig=plt.figure()
            ax = fig.add_subplot(111)
          
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    
                    if minposx == 0:
                        minposx = -1
                    if minposy == 0:
                        minposy = -1
                    if maxposx == 0:
                        maxposx = 1
                    if maxposy == 0:
                        maxposy = 1


                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposy-maxposy/10,maxposy+maxposy/10)
            
            def DefShapeFrame(k,Model,el,scale):
            #Plot deformed shape for a specific time
            #Extrapolate and plot Deformed Shape
                defshape=[]
                for i in range(Model.numel):
                    el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                    Uel=np.zeros((1,Model.nf))
                    for j in range(len(el_boun)):
                        if el_boun[j]==0:
                            Uel[0,j]=Model.Udy[Model.FDOF.index(el[i].id_v[j]-1),k]
                        else:
                            Uel[0,j]=0
                    
                    #This should be a function (Create br)
                    ibr=np.zeros((6,6))
                    r=0
                    for n in range(2):
                        for j in range(len(el[i].ROT)):
                            for p in range(len(el[i].ROT)):
                                ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                        r=r+j+1
                    
                    el_U_l=np.matmul(ibr,scale*np.transpose(Uel))
                    
                    #Start with axial
                    dL=el_U_l[3]-el_U_l[0]
                    
                    #Create Shape Functions
                    xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                    x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[3],20))
                    x=x.squeeze()
                    xl=np.matrix(xl)
                    x=np.matrix(x)
                    l_y=[]
                    
                    for n in range(x.size):
                        N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                        N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                        N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                        N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                        
                        N=np.array([N1,N2,N3,N4])
                        Uy=np.array([el_U_l[1],el_U_l[2],el_U_l[4],el_U_l[5]])
                        
                        l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                        
                        
                    l_y=np.matrix(l_y)
                    
                    def_l=np.concatenate((x,l_y))
                    
                    def_g=np.matmul(el[i].ROT[0:2,0:2],def_l)
                    def_g=def_g.tolist()
                    
                    def_g[0]=[(x+el[i].node[0][0]) for x in def_g[0]]
                    def_g[1]=[(x+el[i].node[0][1]) for x in def_g[1]]
                    
                    defshape.extend(ax.plot(def_g[0],def_g[1],color='r',linewidth=.75))
                
                return defshape
                
            animation = FuncAnimation(fig,func=DefShapeFrame,frames=range(Model.Udy.shape[1]),fargs=(Model,el,scale,),interval=10,blit=True)
            plt.show()
        else:
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
          
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            maxposx=0
            minposx=0
            maxposy=0
            minposy=0
            maxposz=0
            minposz=0
            for i in range(len(el)):
                for j in range(len(el[i].node)):
                    if el[i].node[j][0]>maxposx:
                        maxposx=el[i].node[j][0]
                    if el[i].node[j][0]<minposx:
                        minposx=el[i].node[j][0]
                    if el[i].node[j][1]>maxposy:
                        maxposy=el[i].node[j][1]
                    if el[i].node[j][1]<minposy:
                        minposy=el[i].node[j][1]
                    if el[i].node[j][2]>maxposz:
                        maxposz=el[i].node[j][2]
                    if el[i].node[j][2]<minposz:
                        minposz=el[i].node[j][2]
                            
            ax.set_xlim(minposx-maxposx/10,maxposx+maxposx/10)
            ax.set_ylim(minposz-maxposz/10,maxposz+maxposz/10)
            ax.set_zlim(minposy-maxposy/10,maxposy+maxposy/10)
            
            def DefShapeFrame(k,Model,el,scale):
            #Plot deformed shape for a specific time
            #Extrapolate and plot Deformed Shape
            
                defshape=[]
                for i in range(Model.numel):
                    el_boun=el[i].BOUN[0]+el[i].BOUN[1]
                    Uel=np.zeros((1,Model.nf))
                    for j in range(len(el_boun)):
                        if el_boun[j]==0:
                            Uel[0,j]=Model.Udy[Model.FDOF.index(el[i].id_v[j]-1),k]
                        else:
                            Uel[0,j]=0
                    
                    #This should be a function (Create br)
                    ibr=np.zeros((12,12))
                    r=0
                    for n in range(4):
                        for j in range(len(el[i].ROT)):
                            for p in range(len(el[i].ROT)):
                                ibr[r+p,r+j]=np.linalg.inv(el[i].ROT)[p,j]
                        r=r+j+1
                    
                    el_U_l=np.matmul(ibr,scale*np.transpose(Uel))
                    
                    #Start with axial
                    dL=el_U_l[6]-el_U_l[0]
                    
                    #Create Shape Functions
                    xl=np.transpose(np.linspace(0,el[i].L+dL,20))
                    x=np.transpose(np.linspace(0+el_U_l[0],el[i].L+el_U_l[6],20))
                    x=x.squeeze()
                    x=np.matrix(x)
                    xl=np.matrix(xl)
                    l_z=[]
                    l_y=[]
                    
                    for n in range(x.size):
                        N1=1-3*np.power(xl[0,n],2)/(np.power(el[i].L+dL,2))+(2*np.power(xl[0,n],3)/(np.power(el[i].L+dL,3)))
                        N2=xl[0,n]-(2*np.power(xl[0,n],2)/(el[i].L+dL))+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)  
                        N3=3*np.power(xl[0,n],2)/np.power(el[i].L+dL,2)-2*np.power(xl[0,n],3)/np.power(el[i].L+dL,3)
                        N4=-np.power(xl[0,n],2)/(el[i].L+dL)+np.power(xl[0,n],3)/np.power(el[i].L+dL,2)
                        
                        N=np.array([N1,N2,N3,N4])
                        Uz=np.array([el_U_l[2],el_U_l[4],el_U_l[8],el_U_l[10]])
                        Uy=np.array([el_U_l[1],el_U_l[5],el_U_l[7],el_U_l[11]])
                        
                        l_z.append(np.matmul(np.transpose(N),Uz)[0,0])
                        l_y.append(np.matmul(np.transpose(N),Uy)[0,0])
                        
                        
                    l_z=np.matrix(l_z)
                    l_y=np.matrix(l_y)
                    
                    def_l=np.concatenate((x,l_y,l_z))
                    
                    def_g=np.matmul(el[i].ROT,def_l)
                    def_g=def_g.tolist()
                    
                    def_g[0]=[(x+el[i].node[0][0]) for x in def_g[0]]
                    def_g[1]=[(x+el[i].node[0][1]) for x in def_g[1]]
                    def_g[2]=[(x+el[i].node[0][2]) for x in def_g[2]]
                    
                    defshape.extend(ax.plot(def_g[0],def_g[2],def_g[1],color='r',linewidth=.75))
                
                return defshape
              
            animation = FuncAnimation(fig,func=DefShapeFrame,frames=range(Model.Udy.shape[1]),fargs=(Model,el,scale,),interval=10,blit=True)
            plt.show()
        return animation
    
    def Shear_y(self,Model,el):
        if Model.ModelSet[0]=='Planar':
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])

                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')
            
            
            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')
 
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r',linewidth=.75)
            
            #find shear values and add wl/2 to the shear values
            Model.sheary=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*6)-6):((i+1)*6),0]
                
                #Review Filippou to figure out what this sign convention is about...
                elqy=[elq[1,0],-elq[4,0]]
                
                if el[i].w_y!=0:
                    reacy=[el[i].w_y*el[i].L/2,-el[i].w_y*el[i].L/2]
                else:
                    reacy=[0,0]
                
                elqy=np.array(elqy)-np.array(reacy)
                
                
                x=np.matrix([0,0,el[i].L,el[i].L])
                y=np.matrix([0,maxL/20*elqy[0]/max(abs(np.array(elqy))),maxL/20*elqy[1]/max(abs(np.array(elqy))),0])
        
                coord_l=np.concatenate((x,y))
                coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                
                node_i=np.array(el[i].node[0])
                
                coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[1], color='b',linewidth=.75)
                Model.sheary.append(elqy)
                
                ax.text(coord_g[0][1]\
                        ,coord_g[1][1],str('{:.3f}'.format(elqy[0])))
                ax.text(coord_g[0][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elqy[1])))
        else:
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            z=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                z.append(el[i].node[0][2])
                z.append(el[i].node[1][2])
                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            
            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            BCz=np.zeros(2)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            BCz=np.zeros(10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==2:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.zeros(10)+Model.NODE[i][1]
                            BCz=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=np.zeros(72)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=maxL/40*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==4:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=np.zeros(72)+Model.NODE[i][1]
                            BCz=maxL/40*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
            
            #find shear values and add wl/2 to the shear values
            Model.sheary=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*12)-12):((i+1)*12),0]
                
                #Review Filippou to figure out what this sign convention is about...
                elqy=[elq[1,0],-elq[7,0]]
                
                if el[i].w_y!=0:
                    reacy=[el[i].w_y*el[i].L/2,-el[i].w_y*el[i].L/2]
                else:
                    reacy=[0,0]
                
                elqy=np.array(elqy)-np.array(reacy)
                
                
                x=np.matrix([0,0,el[i].L,el[i].L])
                y=np.matrix([0,maxL/20*elqy[0]/max(abs(np.array(elqy))),maxL/20*elqy[1]/max(abs(np.array(elqy))),0])
                z=np.matrix([0,0,0,0])
        
                coord_l=np.concatenate((x,y,z))
                coord_g=np.matmul(el[i].ROT,coord_l)
                
                node_i=np.array(el[i].node[0])
                
                coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                Model.sheary.append(elqy)
                
                ax.text(coord_g[0][1],coord_g[2][1]\
                        ,coord_g[1][1],str('{:.3f}'.format(elqy[0])))
                ax.text(coord_g[0][-2],coord_g[2][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elqy[1])))
        
    def Shear_z(self,Model,el):
        if Model.ModelSet[0]=='Planar':
            print('ERROR: No z shear in a Planar Analysis')
        else:
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            z=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                z.append(el[i].node[0][2])
                z.append(el[i].node[1][2])
                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            
            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            BCz=np.zeros(2)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            BCz=np.zeros(10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==2:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.zeros(10)+Model.NODE[i][1]
                            BCz=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=np.zeros(72)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=maxL/40*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==4:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=np.zeros(72)+Model.NODE[i][1]
                            BCz=maxL/40*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
            
            #find shear values and add wl/2 to the shear values
            Model.shearz=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*12)-12):((i+1)*12),0]
                
                #Review Filippou to figure out what this sign convention is about...
                elqz=[-elq[2,0],elq[8,0]]
                
                if el[i].w_y!=0:
                    reacz=[-el[i].w_z*el[i].L/2,el[i].w_z*el[i].L/2]
                else:
                    reacz=[0,0]
            
                elqz=np.array(elqz)-np.array(reacz)
                
                
                x=np.matrix([0,0,el[i].L,el[i].L])
                y=np.matrix([0,0,0,0])
                z=np.matrix([0,maxL/20*elqz[0]/max(abs(np.array(elqz))),maxL/20*elqz[1]/max(abs(np.array(elqz))),0])
        
                coord_l=np.concatenate((x,y,z))
                coord_g=np.matmul(el[i].ROT,coord_l)
                
                node_i=np.array(el[i].node[0])
                
                coord_g[:,0]=coord_g[:,0]+np.c_[node_i]
                coord_g[:,1]=coord_g[:,1]+np.c_[node_i]
                coord_g[:,2]=coord_g[:,2]+np.c_[node_i]
                coord_g[:,3]=coord_g[:,3]+np.c_[node_i]
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                Model.shearz.append(elqz)
                
                ax.text(coord_g[0][1],coord_g[2][1]\
                        ,coord_g[1][1],str('{:.3f}'.format(elqz[0])))
                ax.text(coord_g[0][-2],coord_g[2][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elqz[1])))
            
            
    def Moment_y(self,Model,el):
        if Model.ModelSet[0]=='Planar':
            print('ERROR: No y Moment in a Planar Analysis')
        else:
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            z=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                z.append(el[i].node[0][2])
                z.append(el[i].node[1][2])
                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            
            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            BCz=np.zeros(2)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            BCz=np.zeros(10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==2:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.zeros(10)+Model.NODE[i][1]
                            BCz=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=np.zeros(72)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=maxL/40*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==4:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=np.zeros(72)+Model.NODE[i][1]
                            BCz=maxL/40*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
            
            #find moment values and add the simply supported moment to it
            Model.momenty=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*12)-12):((i+1)*12),0]
                
                elmy=[-elq[4,0],elq[10,0]]
                
                if el[i].w_z!=0:
                    x=np.linspace(0,el[i].L,100)
                    ssmoment=el[i].w_z*x/2*(el[i].L-x)
                else:
                    ssmoment=np.zeros(100)
                
                momenty=np.linspace(elmy[0],elmy[1],100)-ssmoment
                momentyp=el[i].L/20*momenty/max(abs(momenty))
                
                
                x=np.matrix(np.hstack((0,np.linspace(0,el[i].L,100),el[i].L)))
                y=np.matrix(np.hstack((0,np.zeros(100),0)))
                z=np.matrix(np.hstack((0,momentyp,0)))
        
                coord_l=np.concatenate((x,y,z))
                coord_g=np.matmul(el[i].ROT,coord_l)
                
                node_i=np.array(el[i].node[0])
                
                for j in range(coord_g.shape[1]):
                    coord_g[:,j]=coord_g[:,j]+np.c_[node_i]
    
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                Model.momenty.append(momenty)
                
                ax.text(coord_g[0][1],coord_g[2][1]\
                        ,coord_g[1][1],str('{:.3f}'.format(elmy[0])))
                ax.text(coord_g[0][-2],coord_g[2][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elmy[1])))
                    
                maxm=0
                for j in range(len(momenty)):
                    if (abs(momenty[j]))>maxm:
                        index=j
                        maxm=abs(momenty[j])
                        maxs=momenty[j]
                
                if index!=0 and index!=99:
                    ax.text(coord_g[0][index],coord_g[2][index]\
                            ,coord_g[1][index],str('{:.3f}'.format(maxs)))
        
    def Moment_z(self,Model,el):
        if Model.ModelSet[0]=='Planar':
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])

                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('y axis')

            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r')

                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            
                            ax.plot(BCx,BCy, color='r',linewidth=.75)
            
            #find moment values and add the simply supported moment to it
            Model.momentz=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*6)-6):((i+1)*6),0]
                
                elmz=[-elq[2,0],elq[5,0]]
                
                if el[i].w_y!=0:
                    x=np.linspace(0,el[i].L,100)
                    ssmoment=el[i].w_y*x/2*(el[i].L-x)
                else:
                    ssmoment=np.zeros(100)
                
                momentz=np.linspace(elmz[0],elmz[1],100)-ssmoment
                momentzp=el[i].L/20*momentz/max(abs(momentz))
                
                
                x=np.matrix(np.hstack((0,np.linspace(0,el[i].L,100),el[i].L)))
                y=np.matrix(np.hstack((0,momentzp,0)))
                
        
                coord_l=np.concatenate((x,y))
                coord_g=np.matmul(el[i].ROT[0:2,0:2],coord_l)
                
                node_i=np.array(el[i].node[0])
                
                for j in range(coord_g.shape[1]):
                    coord_g[:,j]=coord_g[:,j]+np.c_[node_i]
    
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[1], color='b',linewidth=.75)
                Model.momentz.append(momentz)
                
                ax.text(coord_g[0][1]\
                        ,coord_g[1][2],str('{:.3f}'.format(elmz[0])))
                ax.text(coord_g[0][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elmz[1])))
                maxm=0
                for j in range(len(momentz)):
                    if (abs(momentz[j]))>maxm:
                        index=j
                        maxm=abs(momentz[j])
                        maxs=momentz[j]
                
                if index!=0 and index!=99:
                    ax.text(coord_g[0][index]\
                            ,coord_g[1][index],str('{:.3f}'.format(maxs)))
        else:
            #Plotting
            #Plot Model
            
            x=[]
            y=[]
            z=[]
            maxL=0
            for i in range(Model.numel):
                x.append(el[i].node[0][0])
                x.append(el[i].node[1][0])
                y.append(el[i].node[0][1])
                y.append(el[i].node[1][1])
                z.append(el[i].node[0][2])
                z.append(el[i].node[1][2])
                #find maxL element length for scaling
                if el[i].L>maxL:
                    maxL=el[i].L
            
            fig=plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(x,z,y)
            
            ax.set_xlabel('x axis')
            ax.set_ylabel('z axis')
            ax.set_zlabel('y axis')
            
            
            plt.show()
            
            for i in range(len(Model.NODE)):
                for j in range(len(Model.BOUN[0])):
                    if Model.BOUN[i][j]==1:
                        if j==0:
                            BCx=np.linspace(-maxL/40,maxL/40,2)+Model.NODE[i][0]
                            BCy=np.zeros(2)+Model.NODE[i][1]
                            BCz=np.zeros(2)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==1:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][1]
                            BCz=np.zeros(10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==2:
                            BCx=np.zeros(10)+Model.NODE[i][0]
                            BCy=np.zeros(10)+Model.NODE[i][1]
                            BCz=np.linspace(-maxL/40,maxL/40,10)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r')
                        if j==3:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=np.zeros(72)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=maxL/40*np.cos(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==4:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=np.zeros(72)+Model.NODE[i][1]
                            BCz=maxL/40*np.sin(theta)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
                        if j==5:
                            theta=np.linspace(0,2*np.pi,72)
                            
                            BCx=maxL/40*np.cos(theta)+Model.NODE[i][0]
                            BCy=maxL/40*np.sin(theta)+Model.NODE[i][1]
                            BCz=np.zeros(72)+Model.NODE[i][2]
                            
                            ax.plot(BCx,BCz,BCy, color='r',linewidth=.75)
            
            #find moment values and add the simply supported moment to it
            Model.momentz=[]
            for i in range(len(el)):
                elq=Model.q[(((i+1)*12)-12):((i+1)*12),0]
                
                elmz=[-elq[5,0],elq[11,0]]
                
                if el[i].w_y!=0:
                    x=np.linspace(0,el[i].L,100)
                    ssmoment=el[i].w_y*x/2*(el[i].L-x)
                else:
                    ssmoment=np.zeros(100)
                
                momentz=np.linspace(elmz[0],elmz[1],100)-ssmoment
                momentzp=el[i].L/20*momentz/max(abs(momentz))
                
                
                x=np.matrix(np.hstack((0,np.linspace(0,el[i].L,100),el[i].L)))
                y=np.matrix(np.hstack((0,momentzp,0)))
                z=np.matrix(np.hstack((0,np.zeros(100),0)))
                
        
                coord_l=np.concatenate((x,y,z))
                coord_g=np.matmul(el[i].ROT,coord_l)
                
                node_i=np.array(el[i].node[0])
                
                for j in range(coord_g.shape[1]):
                    coord_g[:,j]=coord_g[:,j]+np.c_[node_i]
    
                
                coord_g=coord_g.tolist()
                ax.plot(coord_g[0],coord_g[2],coord_g[1], color='b',linewidth=.75)
                Model.momentz.append(momentz)
                
                ax.text(coord_g[0][1],coord_g[2][1]\
                        ,coord_g[1][2],str('{:.3f}'.format(elmz[0])))
                ax.text(coord_g[0][-2],coord_g[2][-2]\
                        ,coord_g[1][-2],str('{:.3f}'.format(elmz[1])))
                maxm=0
                for j in range(len(momentz)):
                    if (abs(momentz[j]))>maxm:
                        index=j
                        maxm=abs(momentz[j])
                        maxs=momentz[j]
                
                if index!=0 and index!=99:
                    ax.text(coord_g[0][index],coord_g[2][index]\
                            ,coord_g[1][index],str('{:.3f}'.format(maxs)))
                    