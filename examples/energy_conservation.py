#Demonstration of energy conservation (figs. 5 and 6 in arXiv:2302.01893)
import pipic
from pipic_tools import *
import matplotlib.pyplot as plt
import numpy as np
#===========================SIMULATION INITIALIZATION===========================
Temperature = (2/3)*1e-6*electronMass*lightVelocity**2
Density = 1e+18
DebyeLength = sqrt(Temperature/(4*pi*Density*electronCharge**2))
PlasmaPeriod = sqrt(pi*electronMass/(Density*electronCharge**2))
XMin, XMax = -612.4*DebyeLength, 612.4*DebyeLength
FieldAmplitude = 0.001*4*pi*(XMax - XMin)*electronCharge*Density
nx = 32 # number of cells
iterationsPerPlasmaPeriod = 8
timeStep = PlasmaPeriod/iterationsPerPlasmaPeriod
fig_label = 'ec2' 
#---------------------setting solver and simulation region----------------------
if fig_label == 'boris':
    sim = pipic.init(solver='fourier_boris', nx=nx, XMin=XMin, XMax=XMax)
if fig_label == 'boris_no_div_cleaning':
    sim = pipic.init(solver='fourier_boris', nx=nx, XMin=XMin, XMax=XMax)
    sim.fourierSolverSettings(divergenceCleaning = False)
if fig_label == 'ec':
    sim = pipic.init(solver='ec', nx=nx, XMin=XMin, XMax=XMax)
if fig_label == 'ec2':
    sim = pipic.init(solver='ec2', nx=nx, XMin=XMin, XMax=XMax)
#------------------------------adding electrons---------------------------------
@cfunc(type_addParticles)
def density_callback(r, dataDouble, dataInt):# callback function
    return Density # can be any function of coordinate r[0] 
sim.addParticles(name = 'electron', number = nx*100, \
                 charge = -electronCharge, mass = electronMass, \
                 temperature = Temperature, density = density_callback.address)
#---------------------------setting initial field-------------------------------
shift = pi/nx # a shift to mitigate aliasing 
@cfunc(type_fieldLoop)
def setField_callback(ind, r, E, B, dataDouble, dataInt):
    E[0] = FieldAmplitude*sin(2*pi*r[0]/(XMax - XMin) + shift)
sim.fieldLoop(handler=setField_callback.address)
#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)
#-------------preparing output for electron distrribution f(x, px)--------------
xpx_dist = numpy.zeros((64, 128), dtype=numpy.double)
pxLim = 5*sqrt(Temperature*electronMass)
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
#--------------------preparing output of kinetic energy-------------------------
mc2 = electronMass*lightVelocity**2
m2c2 = (electronMass*lightVelocity)**2
kineticEnergy = numpy.zeros((1, ), dtype=numpy.double) 
@cfunc(type_particleLoop)
def kineticEn_cb(r, p, w, id, dataDouble, dataInt):   
    dataDouble[0] += w[0]*mc2*(sqrt(1+(p[0]**2+p[1]**2+p[2]**2)/m2c2)-1)
def getKineticEn(sim):
    kineticEnergy[0] = 0
    sim.particleLoop(name = 'electron', handler = kineticEn_cb.address, \
                     dataDouble = addressOf(kineticEnergy))
    return kineticEnergy[0]
#--------------------preparing output of field energy-------------------------
factor = ((sim.XMax-sim.XMin)/sim.nx)*((sim.YMax-sim.YMin)/sim.ny)* \
        ((sim.ZMax-sim.ZMin)/sim.nz)/(8*pi)
fieldEnergy = numpy.zeros((1, ), dtype=numpy.double) 
@cfunc(type_fieldLoop)
def fieldEn_cb(ind, r, E, B, dataDouble, dataInt):
    dataDouble[0] += factor*(E[0]**2+E[1]**2+E[2]**2+B[0]**2+B[1]**2+B[2]**2)
def getFieldEnergy(sim):
    fieldEnergy[0] = 0
    sim.fieldLoop(handler=fieldEn_cb.address, dataDouble=addressOf(fieldEnergy))
    return fieldEnergy[0]
#===============================SIMULATION======================================
nPlasmaOscillations = 100 # number of plasma oscillations to be simulated
experimentDuration = nPlasmaOscillations*PlasmaPeriod
nIterations = nPlasmaOscillations*iterationsPerPlasmaPeriod
timeP = numpy.zeros((0, ), dtype=numpy.double)
Efield = numpy.zeros((0, ), dtype=numpy.double)
Etot = numpy.zeros((0, ), dtype=numpy.double)
E0 = 0
for i in range(nPlasmaOscillations*iterationsPerPlasmaPeriod):
    print(i, "/", nPlasmaOscillations*iterationsPerPlasmaPeriod)
    Ek = getKineticEn(sim)
    Ef = getFieldEnergy(sim)
    if i == 0:
        E0 = Ek+Ef  
    timeP = np.append(timeP, i*timeStep/PlasmaPeriod)
    Efield = np.append(Efield, Ef/E0)
    Etot = np.append(Etot, (Ek+Ef)/E0)
    print('relative energy deviation = ', ((Ek+Ef)-E0)/E0)
    if Ef > 2*E0:
        break
    sim.advance(timeStep=timeStep)    
#=====================SAVING SIMULATION RESULT TO FILE==========================
Name = fig_label + '_' + str(iterationsPerPlasmaPeriod)
with open(Name + '.npy', 'wb') as f:
    np.save(f, Etot)
    np.save(f, Efield)
    np.save(f, timeP)