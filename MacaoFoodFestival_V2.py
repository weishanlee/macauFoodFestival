# -*- coding: utf-8 -*-
"""
Created on 
@author: weishan_lee

Optimal Layouts of Macao Food Festival
V2: locations of stands are more than number of stands.
"""
from math import sqrt,exp
import numpy as np
import random as rand
import matplotlib.animation as animation
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os.path
# Function Definitions

# Function to calculate the updated total energy
def energyUpdate(i, j, oldEnergy, randomList, r):
    s = oldEnergy
    lList = len(randomList)

    if abs(i-j)==1:
        if i<j:
            if i == 0: 
                iF = randomList[-1]
            else:
                iF   = randomList[i-1]
    
            iC  = randomList[i]
            iN = randomList[i+1]
    
            jF = randomList[j-1]
            jC = randomList[j]
            if j == lList-1: 
                jN = randomList[0]
            else:
                jN = randomList[j+1] 
            
            if iC!=jF or iN!=jC: print("WRONG! iC!=jF or iN!=jC! LINE 57")
            s = s - r[iF,iC] - r[iC,iN] - r[jC,jN]
            s = s + r[iF,jC] + r[jC,iC] + r[iC,jN] # change r[jC,iN] to r[jC,iC]
            
        else: # i>j
            if j == 0: 
                jF = randomList[-1]
            else:
                jF   = randomList[j-1]
    
            jC  = randomList[j]
            jN = randomList[j+1]
    
            iF = randomList[i-1]
            iC = randomList[i]
            if i == lList-1: 
                iN = randomList[0]
            else:
                iN = randomList[i+1] 
            
            if jC!=iF or jN!=iC: print("WRONG! jC!=iF or jN!=iC! LINE 62")
            s = s - r[jF,jC] - r[jC,jN] - r[iC,iN]
            s = s + r[jF,iC] + r[iC,jC] + r[jC,iN]
    else:
        
        if i == 0: 
            iF = randomList[-1]
        else:
            iF   = randomList[i-1]
    
        iC  = randomList[i]
        
        if i == lList-1: 
            iN = randomList[0]
        else:
            iN = randomList[i+1]
    
        if j == 0: 
            jF = randomList[-1]
        else:
            jF   = randomList[j-1]
    
        jC  = randomList[j]
        if j == lList-1: 
            jN = randomList[0]
        else:
            jN = randomList[j+1]
            
        s = s - r[iF,iC] - r[iC,iN] - r[jF,jC] - r[jC,jN]
        s = s + r[iF,jC] + r[jC,iN] + r[jF,iC] + r[iC,jN]  
    return s

# Function to calculate the total energy of a particular state
def energyFunction(r, nLocations):
    energy = 0.0
    
    for i in range(nLocations):
        if r[i,3] == -1: continue
        for j in range(nLocations):
            if j<=i or (r[j,3]==-1): continue
            energy += r[i,3] * r[j,3] / sqrt( (r[i,1]-r[j,1])**2 + (r[i,2]-r[j,2])**2 ) # Coulumb 
            
    return energy

# output of the energy (energy vs time steps)
def outPutEnVSTime(tRecord, energyRecord):
    data = {'tRecord': tRecord,'energyRecord':energyRecord}
    df = pd.DataFrame(data)
    df_file = open('./energyVSTime.csv','w',newline='') 
    df.to_csv(df_file, sep=',', encoding='utf-8',index=False)
    df_file.close()

def outPutAtomsCoor(rr):
    ## wirte optimal stands layout in csv file
    
    df = pd.DataFrame(columns = ['standID','X','Y','charge'])
    df_file = open('./optimalCoorOfAtoms.csv','w',newline='') 

    standID = []
    X = []
    Y = []
    charge = []

    for i in range(nLocations):
        standID += [rr[i,0]]
        X += [rr[i,1]]
        Y += [rr[i,2]]
        charge +=[rr[i,3]]

    df['standID'] = standID
    df['X'] = X
    df['Y'] = Y
    df['charge'] = charge

    df.to_csv(df_file, sep=',', encoding='utf-8', index=False) 
    df_file.close()
    
def plotOptimalAtomLocations(rr, nLocations):
    x = []
    y = []
    #n = [] 
    for i in range(nLocations):
        if r[i,3] != -1:
            #n += [int(r[i,0])]
            x += [r[i,1]]
            y += [r[i,2]]

    fig, ax = plt.subplots()
    ax.title.set_text("Optimal Atoms Locations")

    ax.scatter(x,y,c='blue')

    #for i, txt in enumerate(n):
    #    ax.annotate(txt, (x[i], y[i]))

    ax.ticklabel_format(useOffset=False)
    plt.grid(True)
    plt.xlim(-1+spaceXMin, spaceXMax+1)
    plt.ylim(-1+spaceYMin, spaceYMax+1)
    plt.savefig("optimalLocOfAtoms.eps")     

def writeLog(msg):
    with open('log.txt', 'a+') as the_file:
        print(msg)
        the_file.write(msg)

import os, psutil
# If previous log.txt file exists, remove it.
if os.path.exists("./log.txt"):
    os.remove("./log.txt")
        
def cpu_stats():
    pid = os.getpid()
    py = psutil.Process(pid)
    memory_use = py.memory_info()[0] / 2. ** 30
    return 'Memory: ' + str(np.round(memory_use, 2)) + 'GB\t'

"""# Define grids in space
spaceMin = 0.0
spaceMax = 60
interval = 60
val = np.linspace(spaceMin,spaceMax,interval+1)     # Define space dimension
nLocations = len(val) * len(val)

XY = [[i,j] for i in val for j in val] 
# XY[n][0] refers to x coordinate at the nth coordinates  
# XY[n][1] refers to y coordinate at the nth coordinates"""
    
# Read coordinates and ratings, i.e.,  effctive charges, of the stands
coordsRatings = pd.read_csv("./coorRatings.csv")
nLocations = coordsRatings.shape[0] 
nAtoms =  coordsRatings[coordsRatings['standID']!=-1].shape[0] # nAtoms just means nStands 
spaceXMin = min(coordsRatings.X)
spaceXMax = max(coordsRatings.X)
spaceYMin = min(coordsRatings.Y)
spaceYMax = max(coordsRatings.Y)


## normalize data
#from sklearn import preprocessing
#min_max_scaler = preprocessing.MinMaxScaler()

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

coordsRatings['normRating'] = np.append(scaler.fit_transform(coordsRatings[coordsRatings['standID']!=-1].Ratings.values.reshape(-1, 1)),
                                        coordsRatings[coordsRatings['standID']==-1].Ratings.values.reshape(-1, 1)
                                        )

# Set up the initial configuration
randomList = rand.sample(range(0,nLocations),nAtoms)
#chargeList = np.array([1.0]*nAtoms) #np.random.uniform(-1,1,nAtoms)

msg = "nAtoms = {}\n".format(nAtoms)
writeLog(msg)
msg = "nLocations = {}\n".format(nLocations)
writeLog(msg)
# Define r array
r = np.empty([nLocations,4])
# r[n,0] refers to stand ID
# r[n,1] refers to X coordinates
# r[n,2] refers to Y coordinates
# r[n,3] refers to charge

for i in range(nLocations):
    #r[i,0] = i
    r[i,1] = coordsRatings.X[i] #XY[i][0]
    r[i,2] = coordsRatings.Y[i] #XY[i][1]

nAssignedCharges = 0
for i in range(nLocations):
    if i in randomList:
        r[i,0] = coordsRatings.standID[nAssignedCharges]  # stand ID of the atom
        r[i,3] = coordsRatings.normRating[nAssignedCharges] # charge of the atom
        nAssignedCharges += 1
        #print("i = {}, nAssignedCharges = {}".format(i,nAssignedCharges))
    else:
        r[i,0] = -1   # representing no atom
        r[i,3] = -1   # representing no atom

# plot locations of atoms
x = []
y = []
#n = []
for i in range(nLocations):
    if r[i,3] != -1:
        #n += [int(r[i,0])]
        x += [r[i,1]]
        y += [r[i,2]]

fig, ax = plt.subplots()
ax.title.set_text("Initial Locations of Atoms")

ax.scatter(x,y,c='blue')

#for i, txt in enumerate(n):
#    ax.annotate(txt, (x[i], y[i]))

ax.ticklabel_format(useOffset=False)
plt.grid(True)
plt.xlim(-1+spaceXMin, spaceXMax+1)
plt.ylim(-1+spaceYMin, spaceYMax+1)
plt.savefig("initLocsOfAtoms.eps")

#Calculate the initial energy
energy = energyFunction(r,nLocations)
initEnergy = energy
minEnergy = initEnergy
print("Initial energy = {:.5f}\n".format(initEnergy))

## If you need animation?
animationOption = False
## If you need to record energy vs time step?
energyVsTime = True

# Create a figure window
if animationOption == True:
    
    # Add rAnimation as the collection of r for animation
    rAnimation = []
    rAnimation.append(r[r[:,3]!=-1][:,1:3].tolist()[:])
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                         xlim=(-1+spaceXMin, spaceXMax+1), ylim=(-1+spaceYMin, spaceYMax+1))
    ax.grid()

    #line, = ax.plot([], [], 'k-',lw=2)
    dot, = ax.plot([], [], color='blue', marker='o', markersize=4, markeredgecolor='blue', linestyle='')
    timeStep_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    energy_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    
    def init():
        """initialize animation"""
        #line.set_data([], [])
        dot.set_data([],[])
        timeStep_text.set_text('')
        energy_text.set_text('')
        return dot, timeStep_text, energy_text

# Simulated annealing
# Main loop
Tmax = 1.0
Tmin = 1e-2
tau = 1e4
targetEnergy = -300

tRecord = []
energyRecord = []

t0=0 # setting up the beginning of the time "lump"
tRecord += [0]
energyRecord += [energy]

firstInitial = True

while (energy>targetEnergy):
    
    if firstInitial == False: 
        # Set up another initial configuration
        randomList = rand.sample(range(0, nLocations), nAtoms)

        #Define r array
        r = np.empty([nLocations,4])
        for i in range(nLocations):
            #r[i,0] = i
            r[i,1] = coordsRatings.X[i]
            r[i,2] = coordsRatings.Y[i]

        nAssignedCharges = 0
        for i in range(nLocations):
            if i in randomList:
                r[i,0] = coordsRatings.standID[nAssignedCharges]  # stand ID of the atom
                r[i,3] = coordsRatings.normRating[nAssignedCharges] # charge of the atom
                nAssignedCharges += 1
                #print("i = {}, nAssignedCharges = {}".format(i,nAssignedCharges))
            else:
                r[i,0] = -1   # representing no atom
                r[i,3] = -1   # representing no charge

        #Calculate the initial energy
        energy = energyFunction(r,nLocations)
        
        if animationOption == True:
            rAnimation.append(r[r[:,3]!=-1][:,1:3].tolist()[:])

    T = Tmax
    t = 0
    while (T>Tmin):
        if energy<targetEnergy: break
        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Choose two sites to swap and make sure they are distinct
        i,j = rand.randrange(0,nLocations),rand.randrange(0,nLocations)  # change from 1 to 0
        while ( i==j or ( r[i,3]==-1 and r[j,3]==-1 ) ):
            i,j = rand.randrange(0,nLocations),rand.randrange(0,nLocations)  # change from 1 to 0
                
        # Swap them and calculate the change in energy
        oldEnergy = energy
    
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,3],r[j,3] = r[j,3],r[i,3]
        
        energy = energyFunction(r,nLocations)   
        
        """energy = energyUpdate(i,j,oldEnergy,randomList,r)
        randomList[i], randomList[j] = randomList[j], randomList[i]
        energyCheck = energy(randomList, r)
        if abs(energy-energyCheck)>1e-4:
            randomList[i], randomList[j] = randomList[j], randomList[i]
            msg = "energy Error! Line 315.\n" +\
                  "i = {}, j = {}, randomList[i] = {}, randomList[j] = {}\n".format(i,j,randomList[i],randomList[j]) +\
                  "energy = {}, energyCheck = {}".format(energy,energyCheck)
            writeLog(msg)
            sys.exit()"""
        
        deltaEnergy = energy - oldEnergy

        try:
            ans = np.exp(-deltaEnergy/T)
        except OverflowError:
            if -deltaEnergy/T > 0:
                ans = float('inf')
            else:
                ans = 0.0
    
        # If the move is rejected, swap them back again
        if rand.random() > ans:
            
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,3],r[j,3] = r[j,3],r[i,3]
            
            energy = oldEnergy
            if np.abs(energy - energyFunction(r,nLocations))>1e-5:
                print("energy: {}".format(energy))
                print("energy: {}".format(energy()))
                print("Error Line 385")

        if animationOption == True:
            rAnimation.append(r[r[:,3]!=-1][:,1:3].tolist()[:])
            
        if energyVsTime == True:
            if t%3000==0:
                tRecord += [t0+t]
                energyRecord += [energy]
        
        if energy < minEnergy: 
            minEnergy = energy
            outPutEnVSTime(tRecord, energyRecord)
            outPutAtomsCoor(r)
            dt = datetime.now()
            
            msg = str(dt.year) + '/' + str(dt.month)  + '/' + str(dt.day) + ' ' +\
                  str(dt.hour) + ':' + str(dt.minute) + ':' + str(dt.second) +'\t'
            writeLog(msg)
            writeLog(cpu_stats())
            msg = "Delta energy = {:.5f}\t".format(deltaEnergy)
            writeLog(msg)
            msg = "New energy = {:.5f}\n".format(energy)
            writeLog(msg)    
             
        
    t0 = t0 + t # go to next time "lump"
    firstInitial = False
# End of Main Loop
if animationOption == True:
    def animate(i):
        """perform animation step"""
        xx = []
        yy = []
        for j in range(nAtoms):
            xx.append(rAnimation[i][j][0])
            yy.append(rAnimation[i][j][1])
        #line.set_data(xx,yy)
        dot.set_data(xx,yy)
        timeStep_text.set_text('reduced time = %d ' % tRecord[i] )
        energy_text.set_text('reduce energy = %.3f ' % energyRecord[i] )
        return dot, timeStep_text, energy_text
    
    ani = animation.FuncAnimation(fig, animate, frames=tRecord,
                                  interval=20, blit=True, init_func=init) #interval: duration of time in millisecond between two frames 
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=800, metadata=dict(artist='Me'), bitrate=-1) 
    ani.save('animation.mp4', writer=writer)
    plt.show() # modified on 2021 Dec 20

msg = "The initial energy = {:.5f}\n".format(initEnergy)
writeLog(msg)
msg = "The optimal energy = {:.5f}\n".format(energy)
writeLog(msg)

# plot energy vs t
plt.figure()
plt.title("Energy vs Iteration")
ax = plt.gca()
energyVsTime = pd.read_csv( "./energyVSTime.csv") 
plt.plot(energyVsTime.tRecord,energyVsTime.energyRecord,'k-')
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(5) # number of minor intervals per major # inteval
minorLocatorY = AutoMinorLocator(5)
ax.set_xlabel("Iteration",size = 16)
ax.set_ylabel("Energy",size = 16)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
plt.grid(True)
plt.savefig("energyVsIteration.eps")
plt.show()   

energyCheck = energyFunction(r,nLocations)
msg = "The checked optimal energy = {:.5f}".format(energyCheck)
writeLog(msg)
plotOptimalAtomLocations(r, nLocations)