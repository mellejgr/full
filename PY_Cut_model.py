import csv
from PY_INPUT import *
from datetime import datetime
from scipy.spatial import distance
import os    
import _thread
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.stats import kde
import csv
import math

domain=8 ## divisible by 8
Gmesh=2 ### divisible by 2
print('Block-size='+str(int(domain/Gmesh))+'x'+str(int(domain/Gmesh)))
run="YS"

inp='Nodes_S150'

now = datetime.now().time() # time object
print("Start_time =", now)

if geometry=="triangular":
    ch=round(math.sqrt(3)*l/2, 5)
    cv=l/2
    if (nx % 2) != 0:
        nx=nx+2
    ny=int(nx/2)
    if method=="VTS":
        ny=ny*2

    crackcells=int(nx*aw)
    if (crackcells % 2) == 0:
        cracksize=crackcells*ch+ch/2
    else:
        cracksize=crackcells*ch-ch/2

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

if geometry=="hexagonal":
    cv=l/2             
    ch=round(math.sqrt(3)*l/2, 5)
    
    if (nx % 2) != 0:
        nx=nx+1
    ny=int(nx/2)
    crackcells=int(nx*aw)
    cracksize=crackcells*ch*2
    
##########################################################################
################################### FUNCTIONS ############################
##########################################################################

### FIND NODE CLOSEST TO GIVEN COORDINATES ###
def closest_node(node, nodes): #node= node to search for; nodes=list to search in #returns index
    closest_index = distance.cdist([node], nodes).argmin()
    return closest_index
    

##########################################################################
####################### READ NODES & ELEMENTS ############################
##########################################################################

Xmax=nx*ch/2+domain/2*ch
Xmin=nx*ch/2-domain/2*ch
Ymax=domain*cv
Ymin=-domain*cv

with open('Elements.txt', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)
elements=[[int(x) for x  in sublist] for sublist in data]
for els in elements:
    els.pop(0)


with open('Nodes.txt', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)   
nodes=[[float(x) for x  in sublist] for sublist in data] 
for n in range(0,len(nodes)):
    nodes[n][0]=int(nodes[n][0])   
    

with open(inp+'.txt', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)   
defnodes=[[float(x) for x  in sublist] for sublist in data] 
for n in range(0,len(defnodes)):
    defnodes[n][0]=int(defnodes[n][0]) 


# cut elements and nodes outside of domain:  
addNODE=[]
deleteNODE=[]
defNODE=[]
tolerance=l/(mesh*8)
if len(nodes)==len(defnodes):
    print('Nodes files ok')
for i in range(0,len(nodes)):
    if ((nodes[i][1] <= (Xmax+tolerance) and nodes[i][1]>=(Xmin-tolerance)) and ((Ymin-tolerance) <= nodes[i][2] <= (Ymax+tolerance))):
        addNODE.append(nodes[i])
        defNODE.append(defnodes[i])
    else:
        deleteNODE.append(nodes[i][0])      
nodes=addNODE
defnodes=defNODE

elements =[x for x in elements if not set(x).intersection(deleteNODE)] 


############################################################################################
################################### DIVIDE INTO BLOCKS #####################################
############################################################################################

now = datetime.now().time() # time object
print("Blocks_time =", now)

Xdivisor=(Xmax-Xmin)/(Gmesh)
Ydivisor=(Ymax-Ymin)/(Gmesh)

squares=[]
squares_elements=[]
deformed=[]
for x in range(0,int(Gmesh)):
    x1=Xmin+x*Xdivisor-tolerance
    x2=x1+Xdivisor+2*tolerance
    
    deleteNODE=[]
    for i in range(0,len(nodes)):
        if not (nodes[i][1] <= x2 and nodes[i][1]>=x1):
            deleteNODE.append(nodes[i][0])
    
    #remove elements that are not from Xblock
    subelements =[x for x in elements if not set(x).intersection(deleteNODE)]
    
    
    
    for y in range(0,int(Gmesh)):
        y1=Ymin+y*Ydivisor-tolerance
        y2=y1+Ydivisor+2*tolerance
        

        addNODE=[]
        deleteNODE=[]
        maybe=[]
        defNODE=[]
            
        for i in range(0,len(nodes)):
            if ((nodes[i][1] <= x2 and nodes[i][1]>=x1) and (y1 <= nodes[i][2] <=y2 )):
                addNODE.append(nodes[i])
                defNODE.append(defnodes[i])
                if (nodes[i][2]<tolerance and nodes[i][2]>-tolerance):
                    maybe.append(nodes[i])
            else:
                deleteNODE.append(nodes[i][0])
        
        print(maybe)
        def_squares=defNODE            
        #block to remove duplicates at crack   
        seen = set()
        if y1>0-2*tolerance and y1<0:
            for sl in sorted(maybe, reverse=True):
                if sl[1] not in seen:
                    seen.add(sl[1])
                else:
                    deleteNODE.append(sl[0])
                    addNODE.remove(sl)
        elif y2<0+2*tolerance and y2>0:
            for sl in sorted(maybe):
                if sl[1] not in seen:
                    seen.add(sl[1])
                else:
                    deleteNODE.append(sl[0])
                    addNODE.remove(sl)
        
        
        #remove elements that are not from block
        keep_list =[x for x in subelements if not set(x).intersection(deleteNODE)]
        
        
        deformed.append(def_squares)
        squares.append(addNODE)
        squares_elements.append(keep_list)
        
        
############################################################################################
##################################### CREATE BOUNDARY SETS #################################
############################################################################################

now = datetime.now().time() # time object
print("BC_time =", now)

BC_set=[]
number=0
for i in range(0,Gmesh):
    for j in range(0,int(Gmesh)):       
        number=number+1
        from operator import itemgetter
        Xleft=min(x[1] for x in squares[number-1])
        Ybottom=min(x[2] for x in squares[number-1])
        Xright=max(x[1] for x in squares[number-1])
        Ytop=max(x[2] for x in squares[number-1])

              
        left=[]
        right=[]
        top=[]
        bottom=[]
        seen_b = set()
        seen_t = set()
        for i in squares[number-1]:
            if i[1]==Xleft:
                left.append(i[0])
            elif i[1]==Xright:
                right.append(i[0])
            
            if i[2]==Ybottom:
                if i[1] not in sorted(seen_b):
                    seen_b.add(i[1])
                    bottom.append(i)
                    print(seen_b)
            
            elif i[2]==Ytop:
                if i[1] not in sorted(seen_t, reverse=True):
                    seen_t.add(i[1])
                    top.append(i)
                    
        append=[left,right,top,bottom]
        BC_set.append(append)



############################################################################################
######################################## WRITE INPUT #######################################
############################################################################################

now = datetime.now().time() # time object
print("Write_time =", now)
      
for i in range(0,Gmesh*Gmesh):
    for j in (1,2,3,4):
        with open("E"+str(j)+str(j)+"_"+str(i+1)+".inp","w") as file:
            file.write("** PARTS \n")
            file.write("*Part, name=Part-Model-1 \n")
            file.write("*Node\n")
            for k in deformed[i]:
                file.write(str(k[0])+","+str(k[1])+","+str(k[2])+"\n")
            file.write("*Element, type=B22\n")
            c=0
            for k in squares_elements[i]:
                c=c+1
                file.write(str(c)+","+str(k[0])+","+str(k[1])+","+str(k[2])+"\n")
                
            file.write("*Elset, elset=AllElements\n")
                
            for k in range(0,len(squares_elements[i])):
                file.write(str(k+1)+",\n")                 
                
            file.write("** Section: Section-1  Profile: thickness \n")
            file.write("*Beam Section, elset=AllElements, material=Al, temperature=GRADIENTS, section=RECT \n")
            file.write("1., "+ str(t)+" \n")
            file.write("0.,0.,-1. \n")
            file.write("*End Part \n")
            file.write("** ASSEMBLY \n")
            file.write("*Assembly, name=Assembly \n")
            file.write("*Instance, name=Part-1-1, part=Part-Model-1 \n")
            file.write("*End Instance \n")
            file.write("*Nset, nset=left, instance=Part-1-1 \n")
            for k in BC_set[i][0]:
                file.write(str(k)+",\n")
            c=0
            for k in BC_set[i][0]:
                c=c+1
                file.write("*Nset, nset=L"+str(c)+", instance=Part-1-1 \n")
                file.write(str(k)+"\n")
            file.write("*Nset, nset=right, instance=Part-1-1 \n")
            for k in BC_set[i][1]:
                file.write(str(k)+",\n")
            c=0
            for k in BC_set[i][1]:
                c=c+1
                file.write("*Nset, nset=R"+str(c)+", instance=Part-1-1 \n")
                file.write(str(k)+"\n")
            file.write("*Nset, nset=top, instance=Part-1-1 \n")
            
            
                       
            
            for k in sorted(BC_set[i][2], key = lambda x: x[1]):
                file.write(str(k[0])+",\n")
            c=0
            for k in sorted(BC_set[i][3], key = lambda x: x[1]):
                c=c+1
                file.write("*Nset, nset=T"+str(c)+", instance=Part-1-1 \n")
                file.write(str(k[0])+"\n")
                
            file.write("*Nset, nset=bottom, instance=Part-1-1 \n")
            for k in sorted(BC_set[i][3], key = lambda x: x[1]):
                file.write(str(k[0])+",\n")  
            c=0
            for k in sorted(BC_set[i][2], key = lambda x: x[1]):
                c=c+1
                file.write("*Nset, nset=B"+str(c)+", instance=Part-1-1 \n")
                file.write(str(k[0])+"\n")         
            
            c=0
            if j==1:
                for k in range(0,len(BC_set[i][3])):
                    c=c+1
                    file.write("** Constraint: Constraint-"+str(c)+"v\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("B"+str(c)+", 1, 1. \n")
                    file.write("T"+str(c)+", 1, -1. \n")
                    file.write("** Constraint: Constraint-"+str(c)+"u\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("B"+str(c)+", 2, 1. \n")
                    file.write("T"+str(c)+", 2, 1. \n")
                    file.write("** Constraint: Constraint-"+str(c)+"ur\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("B"+str(c)+", 6, 1. \n")
                    file.write("T"+str(c)+", 6, -1. \n")
                    
                    if c!=1:
                        file.write("** Constraint: Constraint-"+str(c)+"v\n")
                        file.write("*Equation \n")
                        file.write("2 \n")
                        file.write("T"+str(c-1)+", 2, 1. \n")
                        file.write("T"+str(c)+", 2, -1. \n")
                    
                    
            elif j==2:
                for k in range(1,len(BC_set[i][3])):
                    c=c+1
                    file.write("** Constraint: Constraint-"+str(c)+"v\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("L"+str(c)+", 1, 1. \n")
                    file.write("R"+str(c)+", 1, 1. \n")
                    file.write("** Constraint: Constraint-"+str(c)+"u\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("L"+str(c)+", 2, 1. \n")
                    file.write("R"+str(c)+", 2, -1. \n")
                    file.write("** Constraint: Constraint-"+str(c)+"ur\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("L"+str(c)+", 6, 1. \n")
                    file.write("R"+str(c)+", 6, -1. \n")
            
                    file.write("** Constraint: Constraint-"+str(c)+"v\n")
                    file.write("*Equation \n")
                    file.write("2 \n")
                    file.write("R"+str(c)+", 1, 1. \n")
                    file.write("R"+str(c+1)+", 1, -1. \n")
            
#            else:
#                for k in range(1,len(BC_set[i][3])):
#                    c=c+1
#                    file.write("** Constraint: Constraint-"+str(c)+"v\n")
#                    file.write("*Equation \n")
#                    file.write("2 \n")
#                    file.write("L"+str(c)+", 1, 1. \n")
#                    file.write("R"+str(c)+", 1, -1. \n")
#                    file.write("** Constraint: Constraint-"+str(c)+"u\n")
#                    file.write("*Equation \n")
#                    file.write("2 \n")
#                    file.write("L"+str(c)+", 2, 1. \n")
#                    file.write("R"+str(c)+", 2, -1. \n")
#                    file.write("** Constraint: Constraint-"+str(c)+"ur\n")
#                    file.write("*Equation \n")
#                    file.write("2 \n")
#                    file.write("L"+str(c)+", 6, 1. \n")
#                    file.write("R"+str(c)+", 6, 1. \n")
            
#                    file.write("** Constraint: Constraint-"+str(c)+"v\n")
#                    file.write("*Equation \n")
#                    file.write("2 \n")
#                    file.write("R"+str(c)+", 1, 1. \n")
#                    file.write("R"+str(c+1)+", 1, -1. \n")
            
            file.write("*End Assembly \n")
    

            #### MATERIALS ####
        
            file.write("** MATERIALS \n")
            file.write("*Material, name=Al \n")
    
            file.write("*Density \n")
            file.write(str(density)+",\n")
            file.write("*Elastic \n")
            file.write(str(E)+", "+str(poisson)+" \n")
    
            file.write("** ---------------------------------------------------------------- \n")
            file.write("** STEP: Step-1 \n")
        
            #### SOLVER ####
    
            file.write("*Step, name=Step-1, nlgeom=NO \n")
            file.write("*Static \n")
            file.write("1., 1., 1e-05, 1. \n")
    
            
            file.write("** BOUNDARY CONDITIONS \n")
            if j==2:
                file.write("** Name: BOTTOM Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("BOTTOM, 2, 2, -0.1\n")
                file.write("** Name: TOP Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("TOP, 2, 2\n")
            elif j==1: 
                file.write("** Name: LEFT Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("LEFT, 1, 1, -0.1\n")
                file.write("** Name: RIGHT Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("RIGHT, 1, 1\n")
            elif j==3:
                file.write("** Name: BOTTOM Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("BOTTOM, 1, 1, -0.1\n")
                file.write("** Name: TOP Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("TOP, 1, 1\n")
                file.write("** Name: BOTTOM2 Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("BOTTOM, 2, 2 \n")
                file.write("** Name: TOP2 Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("TOP, 2, 2\n")
            else:
                file.write("** Name: RIGHT Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("RIGHT, 2, 2, 0.1\n")
                file.write("** Name: LEFT Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("LEFT, 2, 2\n")
                file.write("** Name: RIGHT2 Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("RIGHT, 1, 1 \n")
                file.write("** Name: LEFT2 Type: Displacement/Rotation \n")
                file.write("*Boundary \n")
                file.write("LEFT, 1, 1\n")

            file.write("** OUTPUT REQUESTS \n")
            file.write("*Restart, write, frequency=0 \n")
            file.write("** FIELD OUTPUT: F-Output-1 \n")
            file.write("*Output, field \n")
            file.write("*Node Output \n")
            file.write("COORD, RF, U, UR \n")
            file.write("*Element Output, directions=YES \n")
            file.write("E, S \n")
            file.write("** HISTORY OUTPUT: H-Output-1 \n")
            file.write("*Output, history, variable=PRESELECT \n")
            file.write("*End Step  \n")

############################################################################################
######################################## SUBMIT JOBS #######################################
############################################################################################


now = datetime.now().time() # time object
print("Submit_time =", now)


if run=="YES":
    
    for i in range(0,Gmesh**2,2):
        print('E_'+str(i+1))            
        os.system('abaqus job=E11_'+str(i+1))
        os.system('abaqus job=E22_'+str(i+1))
        os.system('abaqus job=E33_'+str(i+1))
        os.system('abaqus job=E33_'+str(i+2)+' interactive')
        os.system('abaqus job=E44_'+str(i+1))
        os.system('abaqus job=E44_'+str(i+2))
        os.system('abaqus job=E11_'+str(i+2))
        os.system('abaqus job=E22_'+str(i+2)+' interactive')

            
############################################################################################
####################################### READ RESULTS #######################################
############################################################################################

now = datetime.now().time() # time object
print("Results_time =", now)

os.system('abaqus python PY_GET_res.py '+str(Gmesh)+' '+str(Xdivisor)+' '+str(Ydivisor))


############################################################################################
####################################### CREATE PLOTS #######################################
############################################################################################

now = datetime.now().time() # time object
print("Figures_time =", now)

data_path = 'E.txt'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    data = np.array(list(reader)).astype(float)

x=data[:,0]/l
y=data[:,1]/l

Xmin=np.amin(data[:,0])
Xmax=np.amax(data[:,0])
Ymin=np.amin(data[:,1])
Ymax=np.amax(data[:,1])
z=data[:, 2]/Elat

divX=(x.max()-x.min())/(Gmesh-1)
divY=(y.max()*2)/Gmesh
x=data[:,0]/l-(x.max()-x.min())/2-divX/2

xi, yi = np.linspace(x.min(), x.max(), Gmesh), np.linspace(y.min(), y.max(), Gmesh)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

fig,(ax1)=plt.subplots(figsize=(5,5))

fig=plt.imshow(zi, vmin=0, vmax=1, origin='lower', extent=[x.min()-divX/2, x.max()+divX/2, y.min()-divY/2, y.max()+divY/2], cmap='Oranges')
fig=plt.plot([0,x.min()-divX/2],[0,0], 'k-', lw=1)
cbar=plt.colorbar()
cbar.set_label(r'$E_{d}/E_{lat}$',rotation=90)

ax1.set_ylabel(r'$y/l_0$')
ax1.set_xlabel(r'$x/l_0$')

plt.savefig('E.pdf')
plt.show()
 
 
 
z=-data[:, 7]/G

divX=(x.max()-x.min())/(Gmesh-1)
divY=(y.max()*2)/Gmesh
x=data[:,0]/l-(x.max()-x.min())/2-divX/2

xi, yi = np.linspace(x.min(), x.max(), Gmesh), np.linspace(y.min(), y.max(), Gmesh)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

fig,(ax1)=plt.subplots(figsize=(5,5))

fig=plt.imshow(zi, vmin=0, vmax=1, origin='lower', extent=[x.min()-divX/2, x.max()+divX/2, y.min()-divY/2, y.max()+divY/2], cmap='Oranges')
fig=plt.plot([0,x.min()-divX/2],[0,0], 'k-', lw=1)
cbar=plt.colorbar()
cbar.set_label(r'$G_{d}/G_{lat}$',rotation=90)

ax1.set_ylabel(r'$y/l_0$')
ax1.set_xlabel(r'$x/l_0$')

plt.savefig('G.pdf')
plt.show()

now = datetime.now().time()
print("End_time =", now)





    