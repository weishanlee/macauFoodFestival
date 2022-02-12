# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 09:47:32 2022

@author: Wei-shan

Electric Field Lines for Macau Food Festival

Reference:
    Visualizing a vector field with Matplotlib
https://scipython.com/blog/visualizing-a-vector-field-with-matplotlib/
https://github.com/tomduck/electrostatics/blob/master/examples/dipole.py
"""
#import sys
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.patches import Circle
import pandas as pd
from matplotlib.ticker import AutoMinorLocator


def E(q, r0, x, y):
    """Return the electric field vector E=(Ex,Ey) due to charge q at r0."""
    den = np.hypot(x-r0[0], y-r0[1])**3
    return q * (x - r0[0]) / den, q * (y - r0[1]) / den

"""# Grid of x, y points
nx, ny = 64, 64
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)

# Create a multipole with nq charges of alternating sign, equally spaced
# on the unit circle.
nq = 3 # for simple pole, 2 for dipole, 4 for quadrapole, 8 for octapole         
       # Original code:   2**int(sys.argv[1])
charges = []
for i in range(nq):
    q = i%2 * 2 - 1
    charges.append((q, (np.cos(2*np.pi*i/nq), np.sin(2*np.pi*i/nq))))"""

# Read data
optimalCoorOfAtoms = pd.read_csv( "./optimalCoorOfAtoms.csv")
optimalCoorOfAtoms = optimalCoorOfAtoms[optimalCoorOfAtoms.charge!=-1]
optimalCoorOfAtoms.reset_index(drop=True, inplace=True)
EffQ = optimalCoorOfAtoms.charge
xCoor = optimalCoorOfAtoms.X
yCoor = optimalCoorOfAtoms.Y
nq = len(EffQ)
charges = []
for i in range(nq):
    q = EffQ[i]
    xx = xCoor[i]
    yy = yCoor[i]
    charges.append((q, (xx, yy)))
    
# Grid of x, y points
nx, ny = 1000, 1000

x = np.linspace(min(xCoor), max(xCoor), nx)
y = np.linspace(min(yCoor), max(yCoor), ny)
X, Y = np.meshgrid(x, y)    
    

# Electric field vector, E=(Ex, Ey), as separate components
Ex, Ey = np.zeros((ny, nx)), np.zeros((ny, nx))
for charge in charges:
    ex, ey = E(*charge, x=X, y=Y)
    Ex += ex
    Ey += ey

# plot electric field lines
plt.figure()
plt.title("Electric Field Vectors")
ax = plt.gca()
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(2) # number of minor intervals per major # inteval
minorLocatorY = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis
# Plot the streamlines with an appropriate colormap and arrow style
color = len(EffQ) * np.log(np.hypot(Ex, Ey)) #2 * np.log(np.hypot(Ex, Ey))
ax.streamplot(x, y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='<-', arrowsize=1.5)

# Add filled circles for the charges themselves
#charge_colors = {True: '#aa0000', False: '#0000aa'}
#for q, pos in charges:
#    ax.add_artist(Circle(pos, 0.05, color=charge_colors[q>0]))
ax.scatter(xCoor,yCoor,c=EffQ,cmap=plt.cm.get_cmap('inferno'))
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-1+min(xCoor), 1+max(xCoor))
ax.set_ylim(-1+min(yCoor), 1+max(yCoor))
ax.set_aspect('equal')
plt.grid(True)
plt.savefig("electricFieldLine.eps")
plt.show()

## Plot energy density
plt.figure()
plt.title("Energy Density")
ax = plt.gca()
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(2) # number of minor intervals per major # inteval
minorLocatorY = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

# calculate energy density
energyDensity = Ex**2+Ey**2

pcm=ax.pcolormesh(X,Y,np.log(energyDensity),shading='auto',cmap=plt.cm.get_cmap('plasma'))
cbar = plt.colorbar(pcm)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Logarithmic energy density',rotation=270)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-1+min(xCoor), 1+max(xCoor))
ax.set_ylim(-1+min(yCoor), 1+max(yCoor))
ax.set_aspect('equal')
plt.grid(True)
plt.savefig("energyDensity.eps")
plt.show()

"""## Plot energy density histogram
plt.figure()
plt.title("Energy Density Histogram")
ax = plt.gca()
plt.minorticks_on()
minorLocatorX = AutoMinorLocator(4) # number of minor intervals per major # inteval
minorLocatorY = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minorLocatorX) # add minor ticks on x axis
ax.yaxis.set_minor_locator(minorLocatorY) # add minor ticks on y axis

plt.hist(energyDensity.flatten())
plt.xscale("log")
plt.yscale("log")
ax.set_xlabel('Energy Density')
ax.set_ylabel('Counts')         
plt.grid(True)
plt.savefig("energyDensityHistogram.eps")
plt.show()"""

## plot coordinates of atoms
"""fig, ax = plt.subplots()
ax.title.set_text("Optimal Atoms Locations")
    
sc = ax.scatter(xCoor,yCoor,c=EffQ,cmap=plt.cm.get_cmap('inferno'))
cbar = plt.colorbar(sc)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.set_ylabel('Normalized Effective Charge $Q\prime$',rotation=270)

ax.ticklabel_format(useOffset=False)
plt.grid(True)
plt.xlim(-1+min(xCoor), 1+max(xCoor))
plt.ylim(-1+min(yCoor), 1+max(yCoor))
plt.savefig("optimalLocOfAtoms.eps")     
plt.show()
###############################%% from github################################
Plots field lines for dipole."""

"""import electrostatics
from electrostatics import PointCharge, ElectricField, Potential, GaussianCircle
from electrostatics import finalize_plot

# pylint: disable=invalid-name

XMIN, XMAX = min(xCoor), max(xCoor)#-40, 40
YMIN, YMAX = min(yCoor), max(yCoor)#-30, 30
ZOOM = 6
XOFFSET = 0

electrostatics.init(XMIN, XMAX, YMIN, YMAX, ZOOM, XOFFSET)

# Set up the charges, electric field, and potential
charges = [ PointCharge( EffQ[i],[xCoor[i],yCoor[i]]) for i in range(len(EffQ)) ]
#charges = [PointCharge(1, [-1, 0]),
#           PointCharge(-1, [1, 0])]
field = ElectricField(charges)
potential = Potential(charges)

# Set up the Gaussian surface
g = GaussianCircle(charges[0].x, 0.1)

# Create the field lines
fieldlines = []
for x in g.fluxpoints(field, 12):
    fieldlines.append(field.line(x))
fieldlines.append(field.line([10, 0]))

# Create the vector grid
x, y = np.meshgrid(np.linspace(XMIN/ZOOM+XOFFSET, XMAX/ZOOM+XOFFSET, 41),
                      np.linspace(YMIN/ZOOM, YMAX/ZOOM, 31))
u, v = np.zeros_like(x), np.zeros_like(y)
n, m = x.shape
for i in range(n):
    for j in range(m):
        if any(np.isclose(electrostatics.norm(charge.x-[x[i, j], y[i, j]]),
                             0) for charge in charges):
            u[i, j] = v[i, j] = None
        else:
            mag = field.magnitude([x[i, j], y[i, j]])**(1/5)
            a = field.angle([x[i, j], y[i, j]])
            u[i, j], v[i, j] = mag*np.cos(a), mag*np.sin(a)

## Plotting ##

# Electric field lines and potential contours
fig = plt.figure(figsize=(6, 4.5))
potential.plot()
field.plot()
for fieldline in fieldlines:
    fieldline.plot()
cmap = plt.cm.get_cmap('plasma')
    charge.plot()
for charge in charges:
    charge.plot()
finalize_plot()
#fig.savefig('dipole-field-lines.pdf', transparent=True)

# Field vectors
fig = plt.figure(figsize=(6, 4.5))
plt.quiver(x, y, u, v, pivot='mid', cmap=cmap, scale=35)
for charge in charges:
finalize_plot()
#fig.savefig('dipole-field-vectors.pdf', transparent=True)

plt.show()"""