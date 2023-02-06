#Demonstration of energy conservation (the result is fig. 3 in arXiv:2302.01893)
import pipic
from pipic_tools import *
import matplotlib.pyplot as plt
import math
import numpy as np
#===========================SIMULATION INITIALIZATION===========================
Temperature = 1e-6*electronMass*lightVelocity**2
Density = 1e+18
DebyeLength = math.sqrt(Temperature/(4*math.pi*Density*electronCharge**2))
PlasmaPeriod = math.sqrt(math.pi*electronMass/(Density*electronCharge**2))
XMin, XMax = -500*DebyeLength, 500*DebyeLength
FieldAmplitude = 0.001*4*math.pi*(XMax - XMin)*electronCharge*Density
nx = 32 # number of cells
iterationsPerPlasmaPeriod = 128
timeStep = PlasmaPeriod/iterationsPerPlasmaPeriod
solver = 'ec' # choosing solver to be used
#solver = 'boris'
#---------------------setting solver and simulation region----------------------
if solver == 'boris':
    sim = pipic.boris(nx=nx, XMin=XMin, XMax=XMax)
if solver == 'ec':
    sim = pipic.ec(nx=nx, XMin=XMin, XMax=XMax)
#------------------------------adding electrons---------------------------------
@cfunc(type_addParticles)
def density_callback(r, dataDouble, dataInt):# callback function
    return Density # can be any function of coordinate r[0] 
sim.addParticles(name = 'electron', number = nx*100, \
                 charge = -electronCharge, mass = electronMass, \
                 temperature = Temperature, density = density_callback.address)
#---------------------------setting initial field-------------------------------
@cfunc(type_fieldLoop)
def setField_callback(ind, r, E, B, dataDouble, dataInt):
    E[0] = FieldAmplitude*math.sin(2*math.pi*r[0]/(XMax - XMin))
sim.fieldLoop(handler=setField_callback.address)

#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)
#-------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = numpy.zeros((64, 128), dtype=numpy.double)
pxLim = 5*math.sqrt(Temperature*electronMass)
inv_dx_dpx = (xpx_dist.shape[1]/(XMax - XMin))*(xpx_dist.shape[0]/(2*pxLim))
@cfunc(type_particleLoop)
def xpx_callback(r, p, w, id, dataDouble, dataInt):   
    ix = int(xpx_dist.shape[1]*(r[0] - XMin)/(XMax - XMin))
    iy = int(xpx_dist.shape[0]*0.5*(1 + p[0]/pxLim))
    data = carray(dataDouble, xpx_dist.shape, dtype=numpy.double)
    if iy >= 0 and iy < xpx_dist.shape[0]:
        data[iy, ix] += w[0]*inv_dx_dpx
axs[0].set_title('$\partial N / \partial x \partial p_x$ (s g$^{-1}$cm$^{-2}$)')
axs[0].set(ylabel='$p_x$ (cm g/s)')
axs[0].xaxis.set_ticklabels([])
plot0 = axs[0].imshow(xpx_dist, vmax=3*Density/pxLim, \
                      extent=[XMin,XMax,-pxLim,pxLim], interpolation='none', \
                      aspect='auto', cmap='YlOrBr')
fig.colorbar(plot0, ax=axs[0], location='right')
def plot_xpx():
    xpx_dist.fill(0)
    sim.particleLoop(name = 'electron', handler = xpx_callback.address, \
                     dataDouble = addressOf(xpx_dist))
    plot0.set_data(xpx_dist)
#-------------------------preparing output of Ex(x)-----------------------------
Ex = numpy.zeros((32, ), dtype=numpy.double) 
@cfunc(type_it2r)
def Ex_it2r(it, r, dataDouble, dataInt):
    r[0] = XMin + (it[0] + 0.5)*(XMax - XMin)/Ex.shape[0]
@cfunc(type_field2data)
def get_Ex(it, r, E, B, dataDouble, dataInt):
    dataDouble[it[0]] = E[0]
    
axs[1].set_xlim([XMin, XMax])
axs[1].set(xlabel='$x$ (cm)', ylabel='$E_x$ (cgs units)')
axs[1].set_ylim([-FieldAmplitude, FieldAmplitude])
x_axis = np.linspace(XMin, XMax, Ex.shape[0])
plot_Ex_, = axs[1].plot(x_axis, Ex)
def plot_Ex():
    sim.customFieldLoop(numberOfIterations=Ex.shape[0], it2r=Ex_it2r.address, \
                        field2data=get_Ex.address, dataDouble=addressOf(Ex))
    plot_Ex_.set_ydata(Ex)
#--------------------preparing output of kinetic energy-------------------------
mc2 = electronMass*lightVelocity**2
m2c2 = (electronMass*lightVelocity)**2
kineticEnergy = numpy.zeros((1, ), dtype=numpy.double) 
@cfunc(type_particleLoop)
def kineticEn_cb(r, p, w, id, dataDouble, dataInt):   
    dataDouble[0] += w[0]*mc2*(math.sqrt(1+(p[0]**2+p[1]**2+p[2]**2)/m2c2)-1)
def getKineticEn():
    kineticEnergy[0] = 0
    sim.particleLoop(name = 'electron', handler = kineticEn_cb.address, \
                     dataDouble = addressOf(kineticEnergy))
    return kineticEnergy[0]
#--------------------preparing output of field energy-------------------------
factor = ((sim.XMax-sim.XMin)/sim.nx)*((sim.YMax-sim.YMin)/sim.ny)* \
        ((sim.ZMax-sim.ZMin)/sim.nz)/(8*math.pi)
fieldEnergy = numpy.zeros((1, ), dtype=numpy.double) 
@cfunc(type_fieldLoop)
def fieldEn_cb(ind, r, E, B, dataDouble, dataInt):
    dataDouble[0] += factor*(E[0]**2+E[1]**2+E[2]**2+B[0]**2+B[1]**2+B[2]**2)
def getFieldEnergy():
    fieldEnergy[0] = 0
    sim.fieldLoop(handler=fieldEn_cb.address, dataDouble=addressOf(fieldEnergy))
    return fieldEnergy[0]

#===============================SIMULATION======================================
nPlasmaOscillations = 10 # number of plasma oscillations to be simulated
experimentDuration = nPlasmaOscillations*PlasmaPeriod
nIterations = nPlasmaOscillations*iterationsPerPlasmaPeriod
timeP = numpy.zeros((0, ), dtype=numpy.double)
Efield = numpy.zeros((0, ), dtype=numpy.double)
Etot = numpy.zeros((0, ), dtype=numpy.double)
E0 = 0
for i in range(nPlasmaOscillations*iterationsPerPlasmaPeriod):
    Ek = getKineticEn()
    Ef = getFieldEnergy()
    #print('KinEn =', Ek, 'fieldEn =', Ef, 'total =', Ek+Ef)
    if i == 0:
        E0 = Ek+Ef  
    timeP = np.append(timeP, i*timeStep/PlasmaPeriod)
    Efield = np.append(Efield, Ef/E0)
    Etot = np.append(Etot, (Ek+Ef)/E0)
    #print('relative energy deviation = ', abs(E0-(Ek+Ef))/E0)
    if Ef > 2*E0:
        break
    sim.advance(timeStep=timeStep)

#=====================SAVING SIMULATION RESULT TO FILE==========================
Name = solver + '_' + str(iterationsPerPlasmaPeriod)
with open(Name + '.npy', 'wb') as f:
    np.save(f, Etot)
    np.save(f, Efield)
    np.save(f, timeP)

#==========================PLOTTING SAVED CASES=================================
def addPlot(solver,ippp,ax): #ippp is the number of iterations per plasma period
    name = solver + '_' + str(ippp)
    print('adding ', name)
    with open(name + '.npy', 'rb') as f:
        Etot = np.load(f)
        Efield = np.load(f)
        timeP = np.load(f)    
    cmap = plt.cm.get_cmap('gnuplot')
    col = cmap(0.95*math.log(ippp/2)/math.log(64.0))
    linewidth = 3*math.log(ippp)/math.log(128.0)
    if solver == 'boris':
        linewidth = 0.7
        if ippp == 128:
            linewidth = 3
        cmap = plt.cm.get_cmap('winter')
        col = cmap(1.0*math.log(ippp/8)/math.log(16.0))
    ax.plot(timeP, Etot, color=col, linewidth=linewidth, linestyle=':', \
            marker='.', markersize = 2)
    ax.plot(timeP, Efield, label=solver + ', $T_p/dt=$' + str(ippp), color=col, \
            linewidth=linewidth, linestyle='-', marker='.', markersize = 2)

def plotEn(): # plotting saved cases 
    figEn, axEn = plt.subplots()
    addPlot('boris', 128, axEn) # adding cases to the plot
    addPlot('boris', 16, axEn)
    for i in range(7):
        addPlot('ec', 2**(7-i), axEn)
    axEn.legend()
    axEn.set_ylim([0, 1.25])
    axEn.set_xlim([0, 10])
    axEn.set(ylabel='$E_{field}/E_0$, $E_{tot}/E_0$', xlabel='$t/T_p$')
    figEn.savefig('energy_conservation.png', dpi=300)
#plotEn()


