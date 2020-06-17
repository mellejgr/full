import csv
from PY_INPUT import *
from datetime import datetime
from scipy.spatial import distance
import numpy as np
import matplotlib.pyplot as plt
import os

with open('Manager.txt', 'r') as f:
    lines=f.readlines() 
manager=lines[0].split(',')
stepNR=int(manager[0])

##########################################################################
################################### FUNCTIONS ############################
##########################################################################


##########################################################################
########################### CREATE SDEG FIELD ############################
##########################################################################

now = datetime.now().time() # time object
print("DATA_time =", now)


os.system('abaqus cae nogui=PY_RESULTS_COORD_field.py ')
data_path = 'COORD_field.txt'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    data1 = np.array(list(reader)).astype(float)

os.system('abaqus cae nogui=PY_RESULTS_SDEG_field.py ')

##########################################################################
############################# READ SDEG FIELD ############################
##########################################################################

now = datetime.now().time() # time object
print("PLOT_time =", now)

data_path = 'SDEG_field.txt'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    data = np.array(list(reader)).astype(float)

### save only points with SDEG>0.99
mask=data[:,2]>0.99
data1 = data1[mask]
x1=data1[:,0]
y1=data1[:,1]
#SDEG=data[:,2]

##########################################################################
####################### CREATE A FIT TO CRACK ############################
##########################################################################

with open('Crack_propagation-'str(stepNR+1)+'.txt', 'w') as doc:
    for i in range(len(x1)):
        doc.write(str(stepNR+1)+','+str(x1[i])+','+str(y1[i])+'\n')


data = data[mask]
x=data[:,0]
y=data[:,1]
SDEG=data[:,2]
plt.scatter(x, y)
for i in range(20):
    x=np.append(x,cracksize-ch*i)
    y=np.append(y, 0)
    SDEG=np.append(SDEG, 1)

## create polynomial fit to last 10 ch, and create a tangent line
cutX=max(x)-10*ch
mask=x[:]>cutX

#x1=x1[mask]
#y1=y1[mask]
x=x[mask]
y=y[mask]
SDEG=SDEG[mask]

#print('y=',np.poly1d(np.polyfit(x, y, 2)))
a, b, c = np.polyfit(x, y, 2)
print('polynomial with coeficients:',a,b,c)

#AAA=max(x)
#h=-2*ch
#derivative=(a*(AAA+h)**2+b*(AAA+h)+c-(a*(AAA)**2+b*(AAA)+c))/h
#tan=a*(AAA)**2+b*(AAA)+c+derivative*(x-AAA)

###plot crack orientation

plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)))
#plt.plot(np.unique(x), np.poly1d(np.polyfit(x, tan, 1))(np.unique(x)))
#plt.show()

#### find crack-tip location
cutX_fortip=max(x)-1.5*ch
mask=x[:]>cutX_fortip
x_fortip=x[mask]
y_fortip=y[mask]

X_tip=np.average(x_fortip)
Y_tip=np.average(y_fortip)

a, b = np.polyfit(x, y, 1)
y_crack=(a*x+b)+(a*X_tip+b-Y_tip)
plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y_crack, 1))(np.unique(x)))

xmin=cracksize-ch*20
xmax=cracksize+ch*50
ymin=-cv*50
ymax=+cv*50

axes = plt.gca()
axes.set_xlim([xmin,xmax])
axes.set_ylim([ymin,ymax])


plt.savefig('Crack_'+str(stepNR)+'.pdf')
plt.show()

a, b = np.polyfit(x, y_crack, 1)

angle=math.atan(a)

print(a)
#### write location and angle:
with open('Crack_tip'+'.txt', 'w') as doc:
    doc.write(str(stepNR+1)+','+str(X_tip)+','+str(Y_tip)+','+str(angle)+','+str(a)+','+str(b)+'\n')


print('degrees=',math.degrees(angle))
with open('Manager.txt', 'w') as doc:
    doc.write(str(stepNR+1)+','+str(X_tip)+','+str(Y_tip)+','+str(angle)+','+str(a)+','+str(b)+'\n')



