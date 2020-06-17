import math

pathway='C:/Users/Melle/Desktop/TESTING/Results/' ########UPDATE###########

#### SOLVER #####
SOLVER="EXPLICIT"  #"STANDARD" or "EXPLICIT"
nlgeom="YES"
damage="NO"
method="BLM"
stopper="NO"
recover="NO"
mass_scaling=["NO",10]

#### GEOMETRIC / INPUT ####
aw=0.5
nx=400 #has to be >4, even numbers only (for hexagonal dividing by 2 has to give an uneven number)
mesh=8
Keff=1e7
Mixity=0 #0=mode1, 1=mode2
geometry='triangular'
reldensity=0.3
l=0.1
plane="stress"
B=1 ### stress-biaxiality ratio -> (can be negative); upper boewn limited to when contact occurs contact between crack-surfaces

#### MATERIAL ####
material="BRITTLE" #"BRITTLE" OR "EP"
ES=5e9
poissonS=0.33
density=2700
# GAMMA below t 
Sy=324E6
JC=[Sy, 114E6, 0.42, 1.34] # A(yield strength), B, n, m (temperature), - no strain rate
D=[-0.77, 1.45, 0.47,   0., 0] #D1 initial failure, D2 exponential, -D3 triaxiality, D4 - strain rate, D5 temperature
Tmelt=700
Tglass=650
ref_strain_rate=0.001

#### CONTROLS ####
stopCriterion=["SDEG",0.5]
time=1
write_interval=0.01
velocity=0.01

if SOLVER=="STANDARD":
    damage="NO"

if geometry=='triangular':
    A=2*math.sqrt(3)
if geometry=='hexagonal':
    A=2/math.sqrt(3)
if geometry=='diamond':
    A=2
if geometry=='kagome' or geometry=='kagome2':
    A=math.sqrt(3)

t=l*reldensity/A #thickness direction2 does not influence result == 1
Gamma=Sy**2/(ES*l)#Sy*t/2 #50000 #5*

#### LOAD ####

if Mixity==1:
    K1=0
    K2=Keff
else:
    K1=Keff/(1+math.tan(math.pi*Mixity/2))
    K2=math.tan(math.pi*Mixity/2)*K1

E=ES
poisson=poissonS
if plane=="strain":
    E=ES/(1-poissonS**2)
    poisson=poissonS/(1-poissonS)
GS=ES/(2*(1+poissonS))

if geometry=='triangular':
    Elat=reldensity/3*E
    poissonps=1/3
    if plane=="strain":
        Elat=Elat/(1-poissonS**2/3)
        poissonps=(poissonps*3+poissonS**2)/(3-poissonS**2)
if geometry=='hexagonal':
    Elat=3/2*reldensity**3*E
    poissonps=1
    if plane=="strain":
        Elat=2*Elat/(2-3*poissonS**2*reldensity**2)
        poissonps=(poissonps*2+3*reldensity**2*poissonS**2)/(2-3*reldensity**2*poissonS**2)   
if geometry=='kagome':
    Elat=reldensity/3*E
    poissonps=1/3
    if plane=="strain":
        Elat=Elat/(1-poissonS**2/3)
        poissonps=(poissonps*3+poissonS**2)/(3-poissonS**2)
    
G=Elat/(2*(1+poissonps)) ##check

if plane=="strain":
    kappa=3-4*poissonps
else:
    kappa=(3-poissonps)/(1+poissonps)
    
print(G,kappa, Elat, poissonps)


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
    
Tstress=Keff*B/math.sqrt(math.pi*cracksize)