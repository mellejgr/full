import math
from PY_INPUT import *
from scipy.spatial import distance
from datetime import datetime
import numpy as np
################################### INPUTS ############################

now = datetime.now().time() # time object
print("time =", now)

with open('Manager.txt', 'r') as f:
    lines=f.readlines() 
manager=lines[0].split(',')
stepNR=int(manager[0])
Xtip=float(manager[1])
Ytip=float(manager[2])
gamma_crack=float(manager[3])
a=float(manager[4])
b=float(manager[5])

#stepNR=2

if geometry=="triangular":
    if (crackcells % 2) == 0:
        cracksize=crackcells*ch+ch/2
    else:
        cracksize=crackcells*ch-ch/2
        
    cracktip=[Xtip,Ytip]
    
if geometry=="kagome":
    ch=l/2
    cv=round(l/2*math.sqrt(3), 5)
    
    ny=int(nx/2)
    if method=="VTS":
        ny=ny*2
    if (nx % 2) == 0:
        nx=nx+1
        
    crackcells=int(4*nx*aw)
    cracksize=crackcells*ch-ch/2

M=mesh*2

#cracksize=(nx-1)*ch/2
print(crackcells)
    
##########################################################################
################################### FUNCTIONS ############################
##########################################################################

K1=K1*math.cos(gamma_crack)
K2=K1*math.sin(gamma_crack)

### K-FIELD ###
def displacements(r,theta):
    u1=1/(2*math.sqrt(2*math.pi))*(K1/G)*math.sqrt(r)*(kappa-math.cos(theta))*math.cos(theta/2)+1/(2*math.sqrt(2*math.pi))*(K2/G)*math.sqrt(r)*(kappa+2+math.cos(theta))*math.sin(theta/2)
    u2=1/(2*math.sqrt(2*math.pi))*(K1/G)*math.sqrt(r)*(kappa-math.cos(theta))*math.sin(theta/2)-1/(2*math.sqrt(2*math.pi))*(K2/G)*math.sqrt(r)*(kappa-2+math.cos(theta))*math.cos(theta/2)
    field=[u1,u2]
    #w=-(1+kappa)/(4*math.sqrt(2*math.pi*r)*G)*(K1*math.sin(theta/2)+K2*math.cos(theta/2))
    return field

### angle in coordinates relative to crack tip
def angle_func(XX, YY, crack_angle):
    Xlin=XX*math.cos(crack_angle)+YY*math.sin(crack_angle)
    Ylin=-XX*math.sin(crack_angle)+YY*math.cos(crack_angle)
    alpha=math.atan(Ylin/Xlin)
    
    if (Xlin<0 and Ylin>=0):
        theta=alpha+math.pi
    elif (Xlin<0 and Ylin<0):
        theta=-(math.pi-alpha)
    else:
        theta=alpha
        
    return theta

### angle in coordinates relative to crack tip
def BoundaryCondition(X,Y):
    dat=[X,Y]
    #n=closest_node(dat, NODE)
    n=NODE.index(dat) 
    Y=NODE[n][1]
    X=NODE[n][0]
    #previous
    XX=round(X-previousX, 5)
    YY=round(Y-previousY, 5)
    r=math.sqrt(XX**2+YY**2)                 
    theta=angle_func(XX, YY, previousAngle)

    U=displacements(r,theta)
    deltaU1=U[0]
    deltaU2=U[1]
    
    XX=round(X-cracktip[0], 5)
    YY=round(Y-cracktip[1], 5)

    r=math.sqrt(XX**2+YY**2)                 
    theta=angle_func(XX, YY, gamma_crack)

    U=displacements(r,theta)
    
    U1=(U[0]*(1+time_propagation)-deltaU1)/time_propagation
    U2=(U[1]*(1+time_propagation)-deltaU2)/time_propagation              
                    
    coord.append(dat)
    toADD=[n+1,U1,U2]
    BCset.append(toADD)
    
    Ylin=-XX*math.sin(gamma_crack)+YY*math.cos(gamma_crack)
    if ((cracktip[1]>0 and (Y<=0 and Ylin>=0)) or (cracktip[1]<0 and (Y>=0 and Ylin<=0))):
        BCavoid.append(toADD[0])

### FIND NODE CLOSEST TO GIVEN COORDINATES ###
def closest_node(node, nodes): #node= node to search for; nodes=list to search in #returns index
    closest_index = distance.cdist([node], nodes).argmin()
    return closest_index

### FINDS INTERMEDIATE POINTS BETWEEN P1 AND P2 ###
def intermediates(p1, p2, nb_points): #nb_points is number of points to find
    x_spacing = (p2[0] - p1[0]) / (nb_points + 1)
    y_spacing = (p2[1] - p1[1]) / (nb_points + 1)
    return [[p1[0] + i * x_spacing, p1[1] +  i * y_spacing] 
            for i in range(1, nb_points+1)]

### CHECK WHETHER GIVEN NODE ALREADY EXISTS ###
def check(dat): #dat[x,y]
    if dat[1]!=0 or dat[0]>cracksize: #this avoids concatenating nodes at the crack
        try:
            b=NODE.index(dat) 
        except ValueError:      
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            if add==1: ### append lattice-nodes to list
                lattice_nodes.append(nodeNR)
                if YYY==0 and dat[1]==0:
                    crack_nodes.append(nodeNR)
            return(1)
        else:
            intNR.append(b+1)
            return(0)
    else:
        try:
            if Ycell<0:
                b=bottomNODES.index(dat)
            else:
                b=topNODES.index(dat)
        except ValueError:    
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            if add==1:
                lattice_nodes.append(nodeNR)
                if YYY==0 and dat[1]==0:  ### append crack-nodes to list
                    crack_nodes.append(nodeNR)
            return(1)
        else:
            try:
                b=[i for i, n in enumerate(NODE) if n == dat][1]
            except:
                b=NODE.index(dat)
            intNR.append(b+1)
            return(0)           

##########################################################################
############################# CREATE GEOMETRY ############################
##########################################################################
    
    
if geometry=="triangular":    
    
    ELEMENT=[]
    NODE=[]
    topNODES=[]
    bottomNODES=[]
    lattice_nodes=[]
    crack_nodes=[]
    elementNR=1
    ### create nodes
    nodeNR=1
    
    forN1=[]
    forN4=[]
    forN3=[]
    forN6=[]
    file1 = open("Nodes.txt","w")
    file2 = open("Elements.txt","w")
    for YYY in range(-int(ny/2+1), int(ny/2+1)): ## repeats code for each cell
        Ycell=round(cv*(YYY*2+1), 5)
        prevN1=[]
        prevN2=[]
        forforN4=[]
        forforN3=[]
        forforN6=[]
        for XXX in range(0,int(nx/2)):
            
            ### create lattice-nodes ###
            Xcell=round(ch*(XXX*2+1), 5)
            r=math.sqrt((Xcell-cracksize)**2+Ycell**2)
            if r<50*ch:
                M=mesh*2
            else:
                M=2
            x1=round(Xcell-ch, 5)
            x2=x1
            x3=round(Xcell, 5)
            x4=x3
            x5=x3
            x6=round(Xcell+ch, 5)
            x7=x6
            y1=round(Ycell-cv, 5)
            y6=y1
            y2=round(Ycell+cv, 5)
            y7=y2
            y5=round(Ycell+2*cv, 5)
            if YYY==-1 and Xcell<cracksize-ch:
                y5=round(Ycell+cv, 5)
            y3=round(Ycell-2*cv, 5)
            if YYY==0 and Xcell<cracksize-ch:
                y3=round(Ycell-cv, 5)
            y4=round(Ycell, 5)
            
            ## STRUT 1 ##
            intNR=[]
            add=1
            dat=[x1,y1]
            
            #### if already exists
            if XXX==0:
                if YYY==-int(ny/2+1):
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    nodeNR=nodeNR+1
                else:
                    if YYY==0:
                        NODE.append(dat)
                        if Ycell>0:
                            topNODES.append(dat)
                        else:
                            bottomNODES.append(dat)
                        file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                        intNR.append(nodeNR)
                        nodeNR=nodeNR+1
                    else:
                        intNR.append(forN1[len(forN1)-1])
            else:
                intNR.append(prevN1[len(prevN1)-1])
                
            NODE_1=intNR[len(intNR)-1]
            
            
            
            ### fines nodes on a strut given the mesh-size
            data=intermediates([x1, y1], [x2, y2], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
             
            dat=[x2,y2]
            
            if XXX==0:
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
                
                
                forN1.append(intNR[len(intNR)-1])
            else:
                intNR.append(prevN2[len(prevN2)-1]) 
            
            NODE_2=intNR[len(intNR)-1]

            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")     ### writes elements on the strut##         
            add=1
            
            ## STRUT 2 ##
            
            intNR=[]
            dat=[x1,y1] ### always already exists and equal to previous
            intNR.append(NODE_1)
    
            
            data=intermediates([x1, y1], [x4, y4], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            
            dat=[x4,y4]  ##### if Y>0 always new
            if Ycell>3*cv:
                NODE.append(dat)
                topNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
            else:    
                if YYY==-int(ny/2+1) or Ycell>0:
                    nodeNR=nodeNR+check(dat)
                else:
                    intNR.append(forN4[len(forN4)-1][XXX])
                
            NODE_4=intNR[len(intNR)-1]
            forforN3.append(NODE_4)
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
    
    
            ## STRUT 3 ##
            
            intNR=[]
            dat=[x2,y2]
            intNR.append(NODE_2) ### always already exists and equal to previous
            
            data=intermediates([x2, y2], [x4, y4], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x4,y4] ### always already exists and equal to previous
            intNR.append(NODE_4)
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            ## STRUT 4 ##
            if Ycell>0:
                if YYY!=0 or Xcell<cracksize-ch:
                    intNR=[]
                    dat=[x3,y3]
                    if Xcell<cracksize-ch and YYY==0:
                        NODE.append(dat)
                        if Ycell>0:
                            topNODES.append(dat)
                        else:
                            bottomNODES.append(dat)
                        file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                        intNR.append(nodeNR)
                        nodeNR=nodeNR+1
                    else:
                        intNR.append(forN3[len(forN3)-1][XXX])
                    
                    data=intermediates([x3, y3], [x4, y4], nb_points=M-1)
                    add=0
                    for i in range(0,len(data)):  
                        dat=[round(data[i][0], 6),round(data[i][1], 6)]
                        NODE.append(dat)
                        if Ycell>0:
                            topNODES.append(dat)
                        else:
                            bottomNODES.append(dat)
                        file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                        intNR.append(nodeNR)
                        nodeNR=nodeNR+1
                    add=1
            
                    dat=[x4,y4] ### always already exists and equal to previous
                    intNR.append(NODE_4)
                    
                    for i in range(0,int(M/2)):
                        els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                        elementNR=elementNR+1
                        ELEMENT.append(els)
                        file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
                        
            ## STRUT 5 ##           
            else:         
                intNR=[]
                dat=[x4,y4] ### always already exists and equal to previous
                intNR.append(NODE_4)
                
                data=intermediates([x4, y4], [x5, y5], nb_points=M-1)
                add=0
                for i in range(0,len(data)):  
                    dat=[round(data[i][0], 6),round(data[i][1], 6)]
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    nodeNR=nodeNR+1
                add=1
                
                dat=[x5,y5] ### new if Y<0
                if dat[1]!=0 or dat[0]>cracksize:
                    NODE.append(dat)
        
                    bottomNODES.append(dat)
        
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
        
                    lattice_nodes.append(nodeNR)
                    if YYY==0 and dat[1]==0:
                        crack_nodes.append(nodeNR)
                    nodeNR=nodeNR+1
                else:    
                    NODE.append(dat)
                    bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
    
                    lattice_nodes.append(nodeNR)
                        ### append crack-nodes to list
                    crack_nodes.append(nodeNR)
                    nodeNR=nodeNR+1
                
                forforN4.append(intNR[len(intNR)-1])
                
                for i in range(0,int(M/2)):
                    els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                    elementNR=elementNR+1
                    ELEMENT.append(els)
                    file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
                
            ## STRUT 6 ##
                
            intNR=[]
            dat=[x4,y4] ### always already exists and equal to previous
            intNR.append(NODE_4)
            data=intermediates([x4, y4], [x6, y6], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x6,y6]
            if YYY==-int(ny/2+1) or (YYY==0 and Xcell<cracksize):
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            else:
                intNR.append(forN6[len(forN6)-1][XXX])
            
            
            NODE_6=intNR[len(intNR)-1]
            prevN1.append(NODE_6)
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            ## STRUT 7 ##
                
            intNR=[]
            dat=[x4,y4]  ### always already exists and equal to previous
            intNR.append(NODE_4)
            
            data=intermediates([x4, y4], [x7, y7], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
    
            dat=[x7,y7]  #### always new
            
            
            NODE.append(dat)
            
            
            if Ycell<0:
                bottomNODES.append(dat)
            else:
                topNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
            
            NODE_7=intNR[len(intNR)-1]
            forforN6.append(NODE_7)
            prevN2.append(NODE_7)
            
            
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            if XXX==int(nx/2)-1:
                ## STRUT 8 ##
                    
                intNR=[]
                dat=[x6,y6] ### always already exists and equal to previous
                intNR.append(NODE_6)
                
                data=intermediates([x6, y6], [x7, y7], nb_points=M-1)
                add=0
                for i in range(0,len(data)):  
                    dat=[round(data[i][0], 6),round(data[i][1], 6)]
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    nodeNR=nodeNR+1
                
                dat=[x7,y7]  ### always already exists and equal to previous
                intNR.append(NODE_7)
                
                for i in range(0,int(M/2)):
                    els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                    elementNR=elementNR+1
                    ELEMENT.append(els)
                    file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
        forN4.append(forforN4)
        forN3.append(forforN3)
        forN6.append(forforN6)        
    file1.close()
    file2.close()
    

if geometry=="kagome":
    
    ELEMENT=[]
    NODE=[]
    topNODES=[]
    bottomNODES=[]
    lattice_nodes=[]
    crack_nodes=[]
    elementNR=1
    ### create nodes
    nodeNR=1
    
    file1 = open("Nodes.txt","w")
    file2 = open("Elements.txt","w")
    for YYY in range(-int(ny/2+1), int(ny/2+1)): ## repeats code for each cell
        Ycell=round(cv*(YYY*4+2), 5)
        prevN1=[]
        prevN4=[]
        prevN5=[]
        prevN6=[]
        for XXX in range(0,nx):
            
            ### create lattice-nodes ###
            Xcell=round(4*ch*(XXX)+2*ch, 5)
            
            r=math.sqrt((Xcell-cracksize)**2+Ycell**2)
            if r<50*ch*4:
                M=6
            else:
                M=mesh*2
            
            x1=round(Xcell-2*ch, 5)
            x2=x1
            x3=x1
            x4=x1
            x5=round(Xcell-ch, 5)
            x6=x5
            x7=Xcell
            x8=round(Xcell+ch, 5)
            x9=x8
            x10=round(Xcell+2*ch, 5)
            x11=x10
            if XXX==nx-1:
                x12=x10
                x13=x10
            else:
                x12=round(Xcell+3*ch, 5)
                x13=x12
            
            y1=round(Ycell-2*cv, 5)
            y10=y1
            y2=round(Ycell-cv, 5)
            y5=y2
            y8=y2
            y12=y2
            y7=Ycell
            y3=round(Ycell+cv, 5)
            y6=y3
            y9=y3
            y13=y3
            y4=round(Ycell+2*cv, 5)
            y11=y4
            
            ## STRUT 1 ##
            intNR=[]
            add=1
            dat=[x1,y1]
            
            #### if already exists
            
            if XXX==0:
                if YYY==0 and x1<cracksize:
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    lattice_nodes.append(nodeNR)
                    crack_nodes.append(nodeNR)
                    nodeNR=nodeNR+1
                else:
                    nodeNR=nodeNR+check(dat)   
                    
            else:
                intNR.append(prevN1[len(prevN1)-1])
                
                
                
                
            NODE_1=intNR[len(intNR)-1]
            
            
            
            ### fines nodes on a strut given the mesh-size
            data=intermediates([x1, y1], [x5, y5], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
             
            dat=[x5,y5] ### always neq
            if XXX==0:
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
            else:
                intNR.append(prevN5[len(prevN5)-1])
            NODE_5=intNR[len(intNR)-1]
            
            
            
            
            
            
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")     ### writes elements on the strut##         
            add=1
            
            ## STRUT 2 ##
     
            intNR=[]
            dat=[x5,y5]
            
            intNR.append(NODE_5)
    
            
            data=intermediates([x5, y5], [x7, y7], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            
            dat=[x7,y7]  ##### if Y>0 always new
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
                
            NODE_7=intNR[len(intNR)-1]
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
    
    
            ## STRUT 3 ##
            
            intNR=[]
            dat=[x6,y6]
            if XXX==0:
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
            else:
                intNR.append(prevN6[len(prevN6)-1])
                
            NODE_6=intNR[len(intNR)-1] 
            
            data=intermediates([x6, y6], [x7, y7], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x7,y7] ### always already exists and equal to previous
            intNR.append(NODE_7)
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
           #### STRUT 4
            intNR=[]
            dat=[x4,y4]
            if XXX==0:
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
            else:
                intNR.append(prevN4[len(prevN4)-1])
                
            NODE_4=intNR[len(intNR)-1] 
            data=intermediates([x4, y4], [x6, y6], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x6,y6] ### always already exists and equal to previous
            intNR.append(NODE_6)
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            #### STRUT 5
            intNR=[]
            dat=[x5,y5]
            intNR.append(NODE_5)
            
            data=intermediates([x5, y5], [x8, y8], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x8,y8] ### always already exists and equal to previous
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
            
            NODE_8=intNR[len(intNR)-1] 
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            if XXX==0:
                #### STRUT 6
                intNR=[]
                dat=[x2,y2]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
                    
                NODE_2=intNR[len(intNR)-1] 
                data=intermediates([x2, y2], [x5, y5], nb_points=M-1)
                add=0
                for i in range(0,len(data)):  
                    dat=[round(data[i][0], 6),round(data[i][1], 6)]
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    nodeNR=nodeNR+1
                add=1
                
                dat=[x5,y5] ### always already exists and equal to previous
                intNR.append(NODE_5)
                
                for i in range(0,int(M/2)):
                    els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                    elementNR=elementNR+1
                    ELEMENT.append(els)
                    file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
                
                #### STRUT 6
                intNR=[]
                dat=[x3,y3]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                lattice_nodes.append(nodeNR)
                nodeNR=nodeNR+1
                    
                NODE_3=intNR[len(intNR)-1] 
                
                data=intermediates([x3, y3], [x6, y6], nb_points=M-1)
                add=0
                for i in range(0,len(data)):  
                    dat=[round(data[i][0], 6),round(data[i][1], 6)]
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    nodeNR=nodeNR+1
                add=1
                
                dat=[x6,y6] ### always already exists and equal to previous
                intNR.append(NODE_6)
                
                for i in range(0,int(M/2)):
                    els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                    elementNR=elementNR+1
                    ELEMENT.append(els)
                    file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
                    
            #### STRUT 8
            intNR=[]
            dat=[x6,y6]
            intNR.append(NODE_6)
            
            data=intermediates([x6, y6], [x9, y9], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x9,y9] ### always already exists and equal to previous
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
                
            NODE_9=intNR[len(intNR)-1] 
                
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            
            #### STRUT 9
            intNR=[]
            dat=[x8,y8]
            intNR.append(NODE_8)
            
            data=intermediates([x8, y8], [x7, y7], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x7,y7] ### always already exists and equal to previous
            intNR.append(NODE_7)
                
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            #### STRUT 10
            intNR=[]
            dat=[x9,y9]
            intNR.append(NODE_9)
            
            data=intermediates([x9, y9], [x7, y7], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x7,y7] ### always already exists and equal to previous
            intNR.append(NODE_7)
                
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            #### STRUT 11
            intNR=[]
            dat=[x8,y8]
            intNR.append(NODE_8)
            
            data=intermediates([x8, y8], [x10, y10], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x10,y10] ### always already exists and equal to previous
            
            
            
    
            if YYY==0:
                if x10<cracksize:
                    NODE.append(dat)
                    if Ycell>0:
                        topNODES.append(dat)
                    else:
                        bottomNODES.append(dat)
                    file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                    intNR.append(nodeNR)
                    lattice_nodes.append(nodeNR)
                    crack_nodes.append(nodeNR)
                    nodeNR=nodeNR+1
                else:
                    nodeNR=nodeNR+check(dat) 
            else:
                nodeNR=nodeNR+check(dat) 
            
            NODE_10=intNR[len(intNR)-1] 
                
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            #### STRUT 12
            intNR=[]
            dat=[x8,y8]
            intNR.append(NODE_8)
            
            data=intermediates([x8, y8], [x12, y12], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x12,y12] ### always already exists and equal to previous
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
            
            NODE_12=intNR[len(intNR)-1] 
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            
            #### STRUT 13
            intNR=[]
            dat=[x9,y9]
            intNR.append(NODE_9)
            
            data=intermediates([x9, y9], [x13, y13], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x13,y13] ### always already exists and equal to previous
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
            
            NODE_13=intNR[len(intNR)-1] 
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            #### STRUT 14
            intNR=[]
            dat=[x9,y9]
            intNR.append(NODE_9)
            
            data=intermediates([x9, y9], [x11, y11], nb_points=M-1)
            add=0
            for i in range(0,len(data)):  
                dat=[round(data[i][0], 6),round(data[i][1], 6)]
                NODE.append(dat)
                if Ycell>0:
                    topNODES.append(dat)
                else:
                    bottomNODES.append(dat)
                file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
                intNR.append(nodeNR)
                nodeNR=nodeNR+1
            add=1
            
            dat=[x11,y11] ### always already exists and equal to previous
            NODE.append(dat)
            if Ycell>0:
                topNODES.append(dat)
            else:
                bottomNODES.append(dat)
            file1.write(str(nodeNR)+","+str(dat[0])+","+str(dat[1])+"\n")
            intNR.append(nodeNR)
            lattice_nodes.append(nodeNR)
            nodeNR=nodeNR+1
            
            NODE_11=intNR[len(intNR)-1] 
            
            for i in range(0,int(M/2)):
                els=[elementNR, intNR[i*2],intNR[i*2+1],intNR[i*2+2]]
                elementNR=elementNR+1
                ELEMENT.append(els)
                file2.write(str(els[0])+","+str(els[1])+","+str(els[2])+","+str(els[3])+"\n")
            
            
            prevN1.append(NODE_10)
            prevN5.append(NODE_12)
            prevN6.append(NODE_13)
            prevN4.append(NODE_11)
            
    file1.close()
    file2.close()





############################################################
###################### BC ##################################
############################################################
    
with open('Crack_tip.txt', 'r') as f:
    lines=f.readlines() 
reed=lines[len(lines)-1].split(',')
previousX=float(reed[1])
previousY=float(reed[2])
previousAngle=float(reed[3])

BCavoid=[]
if geometry=="triangular":
        

    if method!="VTS":
        BCset=[]
        coord=[]
        #### left
        X=0
        for YYY in range(-int(ny/2+1), int(ny/2+1)):
            if YYY!=0:
                Y=round(cv*(YYY*2), 5)
                dom=BoundaryCondition(X,Y)
        
        ## right
        X=round(ch*(nx), 5)
        for YYY in range(-int(ny/2+1), int(ny/2+1)):
            Y=round(cv*(YYY*2)+2*cv, 5)
            BoundaryCondition(X,Y)
        
        ###top    
        for XXX in range(0,int(nx/2)):
            Y=round((int(ny/2+1)*2)*cv, 5)
            X=round(XXX*ch*2, 5)
            BoundaryCondition(X,Y)            
            
            Y=round((int(ny/2+1)*2)*cv-cv, 5)
            X=round(X+ch, 5)
            BoundaryCondition(X,Y)
                
        ##bottom
        for XXX in range(0,int(nx/2)):
            Y=round(-(int(ny/2+1)*2)*cv, 5)
            X=round(XXX*ch*2+2*ch, 5)
            BoundaryCondition(X,Y)        
            
            Y=round(-(int(ny/2+1)*2)*cv+cv, 5)
            X=round(XXX*ch*2+ch, 5)
            BoundaryCondition(X,Y)
            
                
    else:
        topNodes=[]
        for XXX in range(0,int(nx/2)):
            Y=round((int(ny/2+1)*2)*cv, 5)
            X=round(XXX*ch*2, 5)
            dat=[X,Y]
            n=NODE.index(dat)+1 
            topNodes.append(n)
            
            Y=round((int(ny/2+1)*2)*cv-cv, 5)
            X=round(X+ch, 5)
            dat=[X,Y]
            n=NODE.index(dat) +1
            topNodes.append(n)
            if XXX==int(nx/2)-1:
                Y=round((int(ny/2+1)*2)*cv, 5)
                X=round(X+ch, 5)
                dat=[X,Y]
                n=NODE.index(dat) +1
                topNodes.append(n)
        bottomNodes=[]
        for XXX in range(0,int(nx/2)):
            Y=-round((int(ny/2+1)*2)*cv, 5)
            X=round(XXX*ch*2, 5)
            dat=[X,Y]
            n=NODE.index(dat)+1
            bottomNodes.append(n)
            
            Y=-round((int(ny/2+1)*2)*cv-cv, 5)
            X=round(X+ch, 5)
            dat=[X,Y]
            n=NODE.index(dat) +1
            bottomNodes.append(n)
            if XXX==int(nx/2)-1:
                Y=-round((int(ny/2+1)*2)*cv, 5)
                X=round(X+ch, 5)
                dat=[X,Y]
                n=NODE.index(dat) +1
                bottomNodes.append(n)
    
if geometry=="kagome":
    
    if method!="VTS":
        BCset=[]
        coord=[]
        #### left
        X=0
        
        for YYY in range(-int(ny/2+1), int(ny/2+1)):
            if YYY!=0:
                Y=round(cv*(YYY*4), 5)
                dat=[X,Y]
                #n=closest_node(dat, NODE)
                n=NODE.index(dat) 
                Y=NODE[n][1]
                X=NODE[n][0]
                XX=round(X-cracksize, 5)
                
                r=math.sqrt(XX**2+Y**2)       
                theta=angle_func(XX, YY)
                
                U=displacements(r,theta)
                coord.append(dat)
                toADD=[n+1,U[0],U[1]]
                BCset.append(toADD)
                
            Y=round(cv*(YYY*4+1), 5)
            dat=[X,Y]
            #n=closest_node(dat, NODE)
            n=NODE.index(dat) 
            Y=NODE[n][1]
            X=NODE[n][0]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)          
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
            
            Y=round(cv*(YYY*4+3), 5)
            dat=[X,Y]
            #n=closest_node(dat, NODE)
            n=NODE.index(dat) 
            Y=NODE[n][1]
            X=NODE[n][0]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)          
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
        
        ## right
        X=round(4*ch*nx, 5)
        for YYY in range(-int(ny/2+1), int(ny/2+1)):
            Y=round(cv*(YYY*4+1), 5)
            dat=[X,Y]
            #n=closest_node(dat, NODE)
            n=NODE.index(dat) 
            Y=NODE[n][1]
            X=NODE[n][0]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)          
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
            
            Y=round(cv*(YYY*4+3), 5)
            dat=[X,Y]
            #n=closest_node(dat, NODE)
            n=NODE.index(dat) 
            Y=NODE[n][1]
            X=NODE[n][0]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)          
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
            
            Y=round(cv*(YYY*4+4), 5)
            dat=[X,Y]
            #n=closest_node(dat, NODE)
            n=NODE.index(dat) 
            Y=NODE[n][1]
            X=NODE[n][0]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)          
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
        
        ###top
        
        
        for XXX in range(0,int(nx)):
            Y=round(cv*(ny/2*4+2), 5)
            X=round(XXX*ch*4, 5)
            dat=[X,Y]
            n=NODE.index(dat) 
            X=NODE[n][0]
            Y=NODE[n][1]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)  
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
            
            
        ##bottom
        
        for XXX in range(0,int(nx-1)):
            Y=round(-cv*(ny/2*4+2), 5)
            X=round(XXX*ch*4+4*ch, 5)
            dat=[X,Y]
            n=NODE.index(dat) 
            X=NODE[n][0]
            Y=NODE[n][1]
            XX=round(X-cracksize, 5)
            
            r=math.sqrt(XX**2+Y**2)
            theta=angle_func(XX, YY)
            
            U=displacements(r,theta)  
            coord.append(dat)
            toADD=[n+1,U[0],U[1]]
            BCset.append(toADD)
    else:
        topNodes=[]
        bottomNodes=[]
        Y=round(cv*(ny/2*4+4), 5)
        for XXX in range(0,int(nx)):
            X=round(XXX*ch*4, 5)
            dat=[X,Y]
            n=NODE.index(dat) +1
            topNodes.append(n)
            if XXX==nx-1:
                X=round(XXX*ch*4+4*ch, 5)
                dat=[X,Y]
                n=NODE.index(dat) +1
                topNodes.append(n)
        Y=-round(cv*(ny/2*4+4), 5)
        for XXX in range(0,int(nx)):
            X=round(XXX*ch*4, 5)
            dat=[X,Y]
            n=NODE.index(dat) +1
            bottomNodes.append(n)
            if XXX==nx-1:
                X=round(XXX*ch*4+4*ch, 5)
                dat=[X,Y]
                n=NODE.index(dat) +1
                bottomNodes.append(n)
            

##### write sets #######
    
with open("Sets.txt","w") as file3:
    file3.write("*Nset, nset=Lattice_Nodes \n")
    for i in range(0,len(lattice_nodes)-1):
        file3.write(str(lattice_nodes[i])+",\n")
    file3.write(str(lattice_nodes[len(lattice_nodes)-1])+"\n")


    file3.write("*Nset, nset=Crack_Nodes \n")
    for i in range(0,len(crack_nodes)-1):
        file3.write(str(crack_nodes[i])+",\n")
    file3.write(str(crack_nodes[len(crack_nodes)-1])+"\n")
    
    if method=="VTS":
        file3.write("*Nset, nset=top_Nodes \n")
        for i in range(0,len(topNodes)-1):
            file3.write(str(topNodes[i])+",\n")
        file3.write(str(topNodes[len(topNodes)-1])+"\n")
        
        file3.write("*Nset, nset=bottom_Nodes \n")
        for i in range(0,len(bottomNodes)-1):
            file3.write(str(bottomNodes[i])+",\n")
        file3.write(str(bottomNodes[len(bottomNodes)-1])+"\n")

    file3.write("*Elset, elset=AllElements \n")
    for i in range(0,len(ELEMENT)-1):
        file3.write(str(i+1)+",\n")
    file3.write(str(len(ELEMENT))+"\n")
    
BCcreator=np.array(BCset) 

#print(BCcreator)
#print(' ')
#
#print(' ')
#
#print(' ')
#print(BCset)
if method=="BLM":
    file5 = open("BCs_sets.txt","w")
    for i in range(0,len(BCset)):
        file5.write("*Nset, nset=BCset"+str(i+1)+", instance=Part-1-1 \n" + str(BCset[i][0])+"\n")
    file5.close()
    
    file4 = open("BCs_step"+str(stepNR)+".txt","w")
    j=1
    for i in BCcreator:
        file4.write("**Name: BC_"+str(j)+" Type: Velocity/Angular velocity \n")
        if i[0] not in BCavoid:
            file4.write("*Boundary, op=NEW, type=VELOCITY \n" + "BCset"+str(j)+", 1, 1, "+str(i[1])
                        +"\n" + "BCset"+str(j)+", 2, 2, "+str(i[2])+"\n")
        else:
            file4.write("*Boundary, op=NEW\n")
        j=j+1
    file4.close()


##############################################################
####################### WRITE INPUT FILE #####################
##############################################################

if geometry=="triangular":
    jobname="T"+str(int(reldensity*1000))+"_S"+str(nx)+"_M"+str(mesh)
if geometry=="kagome":
    jobname="K"+str(int(reldensity*1000))+"_S"+str(nx)+"_M"+str(mesh)

with open(str(jobname)+"_step"+str(stepNR)+".inp","w") as INP:

    INP.write("*Heading \n")
    INP.write("*RESTART, READ, STEP="+str(stepNR-1)+", INTERVAL=1, END STEP \n")  
    
    INP.write("** ---------------------------------------------------------------- \n")
    INP.write("** STEP: Step-"+str(stepNR)+" \n")
    
    #### SOLVER ####
    
    if SOLVER=="STANDARD":
        INP.write("*Step, name=Step-"+str(stepNR)+", nlgeom=YES, inc=100000 \n")
        INP.write("*Static, stabilize=2e-05, allsdtol=0, continue=NO \n")
        INP.write("0.01, 1., 1e-08, 0.01 \n")
    else:
        INP.write("*Step, name=Step-"+str(stepNR)+", nlgeom=YES \n")
        INP.write("*Dynamic, Explicit \n")
        INP.write(", "+str(time_propagation)+",, 0.01 \n")
        INP.write("*Bulk Viscosity \n")
        INP.write("0.06, 1.2 \n")        
       
        
        
    INP.write("** BOUNDARY CONDITIONS \n")

    with open("BCs_step"+str(stepNR)+".txt","r") as file:
        INP.write(file.read())
    if SOLVER=="STANDARD":
        INP.write("** CONTROLS \n")
        INP.write("*Controls, reset \n")
        INP.write("*Controls, parameters=time incrementation \n")
        INP.write(", , , , , , , 50, , ,  \n")
    
    #### OUTPUT REQUESTS ####
        
    INP.write("** OUTPUT REQUESTS \n")
    if SOLVER=="STANDARD":
        INP.write("*Restart, write, frequency=0 \n")
        INP.write("** FIELD OUTPUT: F-Output-"+str(stepNR)+" \n")
        INP.write("*Output, field \n")
    else:
        INP.write("*Restart, write, number interval=1, time marks=NO\n")
        INP.write("** FIELD OUTPUT: F-Output-"+str(stepNR)+" \n")
        INP.write("*Output, field, time interval=" +str(time_propagation)+" \n")
        
    INP.write("*Node Output \n")
    INP.write("COORD, U, UR \n")
    INP.write("*Element Output, directions=YES \n")
    INP.write("E, S, SF, PEEQ, SDEG, DMICRT, TRIAX \n")
#    INP.write("** HISTORY OUTPUT: H-Output-"+str(stepNR)+" \n")
#    INP.write("*Output, history, variable=PRESELECT \n")

    INP.write("*End Step  \n")

now = datetime.now().time() # time object
print("time =", now)
print("finished")