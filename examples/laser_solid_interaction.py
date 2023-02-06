import pipic
from pipic_tools import *
import matplotlib.pyplot as plt
import matplotlib.colors as plt_col
from mpl_toolkits.axes_grid1 import make_axes_locatable
from math import *
import numpy as np
#===========================SIMULATION INITIALIZATION===========================
wavelength = 1e-4
nx, XMin, XMax = 128, -8*wavelength, 8*wavelength
ny, YMin, YMax = 64, -4*wavelength, 4*wavelength
dx, dy = (XMax - XMin)/nx, (YMax - YMin)/ny
timeStep = 0.5*dx/lightVelocity
figStride = int(nx/16)

@cfunc(types.double(types.double,types.double,types.double)) #auxiliary function
def cos2shape(x, plateauSize, transitionSize):
    return (abs(x) < plateauSize/2 + transitionSize)* \
           (1 - (abs(x) > plateauSize/2)* \
           sin(0.5*pi*(abs(x) - plateauSize/2)/transitionSize)**2)

#---------------------setting solver and simulation region----------------------
sim = pipic.ec(nx=nx, ny=ny, XMin=XMin, XMax=XMax, YMin=YMin, YMax=YMax) 
#---------------------------setting field of the pulse--------------------------
amplitude_a0 = 100
omega = 2*pi*lightVelocity/wavelength
a0 = electronMass*lightVelocity*omega/(-electronCharge)
fieldAmplitude = amplitude_a0*a0
arrivalDelay = 6*2*pi/omega
incidenceAngle = pi/3
# setup for a basic field generator/absorber
boundarySize = wavelength
supRate = 1 - exp(log(0.01)/(wavelength/(lightVelocity*timeStep)))

@cfunc(types.double(types.double, types.double)) # longitudinal shape
def pulse(eta, phi):
    return (abs(eta-2*pi)< 2*pi)* \
           sin(eta/4)**3*(cos(eta/4)*sin(eta+phi)+sin(eta/4)*cos(eta+phi))

@cfunc(type_fieldLoop)
def field_callback(ind, r, E, B, dataDouble, dataInt):
    if dataInt[0] == 0 or YMax - r[1] < boundarySize:
        # computing longitudinal and transverse coordinates:
        lCoord = sin(incidenceAngle)*r[0] - cos(incidenceAngle)*r[1] + \
                 + (arrivalDelay - dataInt[0]*timeStep)*lightVelocity
        tCoord = cos(incidenceAngle)*r[0] + sin(incidenceAngle)*r[1] 
        tShape = cos2shape(tCoord, 3*wavelength, 2*wavelength)
        # computing both field components of the CP pulse
        f_1 = fieldAmplitude*tShape*pulse(2*pi*lCoord/wavelength + 4*pi, pi)
        f_2 = fieldAmplitude*tShape*pulse(2*pi*lCoord/wavelength + 4*pi, pi/2)
        Ex = cos(incidenceAngle)*f_1
        Ey = sin(incidenceAngle)*f_1
        Bz = f_1
        Bx = -cos(incidenceAngle)*f_2
        By = -sin(incidenceAngle)*f_2
        Ez = f_2
        if dataInt[0] == 0: # setting initial field
            E[0], E[1], E[2], B[0], B[1], B[2] = Ex, Ey, Ez, Bx, By, Bz
        else: # soft reduction of the field difference (pulse generation)
            rate = supRate*cos2shape(YMax - r[1], 0, boundarySize)
            E[0] = (1 - rate)*E[0] + rate*Ex
            E[1] = (1 - rate)*E[1] + rate*Ey
            E[2] = (1 - rate)*E[2] + rate*Ez
            B[0] = (1 - rate)*B[0] + rate*Bx
            B[1] = (1 - rate)*B[1] + rate*By
            B[2] = (1 - rate)*B[2] + rate*Bz
    if r[1] - YMin < boundarySize: # basic absorber
        rate = supRate*cos2shape(r[1] - YMin, 0, boundarySize)
        E[0] = (1 - rate)*E[0]
        E[1] = (1 - rate)*E[1]
        E[2] = (1 - rate)*E[2]
        B[0] = (1 - rate)*B[0]
        B[1] = (1 - rate)*B[1]
        B[2] = (1 - rate)*B[2]
#-----------------------------seting plasma-------------------------------------
N_cr = electronMass*omega**2/(4*pi*electronCharge**2)
density = 100*N_cr
DebyeLength = 0.1*wavelength/64.0
Temperature = 4*pi*density*(electronCharge**2)*DebyeLength**2
particlesPerCell = 100
@cfunc(type_addParticles)
def density_callback(r, dataDouble, dataInt):# callback function 
    return density*cos2shape(r[1] + wavelength, wavelength, wavelength) 
sim.addParticles(name = 'electron', number = int(nx*ny*0.25*particlesPerCell), \
                 charge = -electronCharge, mass = electronMass, \
                 temperature = Temperature, density = density_callback.address)
#=================================OUTPUT========================================
fig, axs = plt.subplots(2, constrained_layout=True)
ax = axs[0]
#------------------------------field output-------------------------------------
tmp = numpy.zeros((1, 1), dtype=numpy.double) # null plot to show the color bar
G_cmap=plt_col.LinearSegmentedColormap.from_list('N',[(1,1,1),(0,0.5,0)],N=32)
cbar = ax.imshow(tmp, vmin=0, vmax=2, \
                 extent=[XMin,XMax,YMin,YMax], interpolation='none', \
                 aspect='equal', origin='lower', cmap=G_cmap)
oBz = numpy.zeros((ny, nx), dtype=numpy.double)
maxField = 2*fieldAmplitude
plot0 = ax.imshow(oBz, vmin=-maxField, vmax=maxField, \
                      extent=[XMin,XMax,YMin,YMax], interpolation='none', \
                      aspect='equal', cmap='RdBu', origin='lower')
ax.set(xlabel='$x$ (cm)', ylabel='$y$ (cm)')
ax.ticklabel_format(axis='both', scilimits=(0,0), useMathText=True)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.05)
fig.add_axes(cax)
bar0 = fig.colorbar(plot0, cax = cax, orientation = "vertical")
bar0.set_label(label='$B_z$ (CGS)')
@cfunc(type_fieldLoop)
def fieldBz_cb(ind, r, E, B, dataDouble, dataInt):
    Bz = carray(dataDouble, oBz.shape, dtype=numpy.double)
    Bz[ind[1], ind[0]] = B[2]
def plot_Field():
    sim.fieldLoop(handler=fieldBz_cb.address, dataDouble=addressOf(oBz))
    plot0.set_data(oBz)
#----------------------------density output-------------------------------------
oN = numpy.zeros((ny, nx), dtype=numpy.double)
@cfunc(type_particleLoop)
def density_cb(r, p, w, id, dataDouble, dataInt):   
    ix = int(oN.shape[1]*(r[0] - XMin)/(XMax - XMin))
    iy = int(oN.shape[0]*(r[1] - YMin)/(YMax - YMin))
    data = carray(dataDouble, oN.shape, dtype=numpy.double)
    if iy >= 0 and iy < oN.shape[0] and ix >= 0 and ix < oN.shape[1]:
        data[iy, ix] += w[0]/(dx*dy*2*density)
N_cmap=plt_col.LinearSegmentedColormap.from_list('N',[(0,0.5,0),(0,0.5,0)],N=9)
plot1 = ax.imshow(oN, vmin=0, vmax=1, alpha=oN, \
                      extent=[XMin,XMax,YMin,YMax], interpolation='none', \
                      aspect='equal', origin='lower', cmap=N_cmap)
ax.set(xlabel='$x$ (cm)', ylabel='$y$ (cm)')
ax.ticklabel_format(axis='both', scilimits=(0,0), useMathText=True)
cax1 = divider.append_axes("right", size="2%", pad=0.85)
fig.add_axes(cax1)
bar1 = fig.colorbar(cbar, cax = cax1, orientation = "vertical")
bar1.set_label(label='$N/N_0$')
def plot_Density():
    oN.fill(0)
    sim.particleLoop(name = 'electron', handler = density_cb.address, \
                     dataDouble = addressOf(oN))
    plot1.set_data(oN)
#----------------------------1dfield output-------------------------------------
oEp = numpy.zeros((2*nx,),dtype=numpy.double) #field P-component of the pulse
oEs = numpy.zeros((2*nx,),dtype=numpy.double) #field S-component of the pulse
EpSize = 8*wavelength # size of the output
@cfunc(type_it2r)
def E_it2r(it, r, dataDouble, dataInt):
    lCoord = (it[0] + 0.5)*EpSize/oEp.shape[0]
    r[0] = lCoord*sin(incidenceAngle)
    r[1] = lCoord*cos(incidenceAngle)
    r[2] = 0
@cfunc(type_field2data)
def get_Ep(it, r, E, B, dataDouble, dataInt):
    dataDouble[it[0]] = -E[0]*cos(incidenceAngle) + E[1]*sin(incidenceAngle)
@cfunc(type_field2data)
def get_Es(it, r, E, B, dataDouble, dataInt):
    dataDouble[it[0]] = E[2]
axs[1].set_xlim([0, EpSize])
axs[1].set_ylim([-maxField, maxField])
axs[1].set(xlabel='$x$ (cm)', ylabel='$E_p$ (cgs units)')
x_axis = np.linspace(0, EpSize, oEp.shape[0])
sim.customFieldLoop(oEp.shape[0],E_it2r.address,get_Ep.address,addressOf(oEp))
sim.customFieldLoop(oEp.shape[0],E_it2r.address,get_Ep.address,addressOf(oEs))
plot_Ep, = axs[1].plot(x_axis, oEp)
#plot_Es, = axs[1].plot(x_axis, oEs)
with open('res_.npy', 'rb') as f: # reading data of RES computation
    res_X_0 = np.load(f)
    res_Ep = -np.load(f)
    res_Es = -np.load(f)
res_X = res_X_0 + -arrivalDelay*lightVelocity
plot_res_Ep, = axs[1].plot(res_X, res_Ep)
#plot_res_Es, = axs[1].plot(res_X, res_Es)

def plot_E(i):
    sim.customFieldLoop(numberOfIterations=oEp.shape[0], it2r=E_it2r.address, \
                        field2data=get_Ep.address, dataDouble=addressOf(oEp))
    sim.customFieldLoop(numberOfIterations=oEs.shape[0], it2r=E_it2r.address, \
                        field2data=get_Es.address, dataDouble=addressOf(oEs))
    plot_Ep.set_ydata(oEp)
    #plot_Es.set_ydata(oEs)
    res_X = res_X_0 + ((-arrivalDelay + i*timeStep)*lightVelocity + \
            wavelength*cos(incidenceAngle)) # due to surface being at 0.5 \mu m
    plot_res_Ep.set_xdata(res_X)
    if i == 21*figStride:
        with open('im21_Ep_' + str(nx) + '.npy', 'wb') as f:
            np.save(f, x_axis)
            np.save(f, oEp)
            #np.save(f, Ez_out)
    #plot_res_Es.set_xdata(res_X)

#===============================SIMULATION======================================
dataInt = numpy.zeros((1, ), dtype=int) # data for passing the iteration number
for i in range(figStride*21 + 1):
    print(i, '/', figStride*21 + 1)
    dataInt[0] = i
    sim.fieldLoop(handler=field_callback.address, dataInt=addressOf(dataInt), \
                  useOmp=True)
    if i%figStride == 0:
        plot_Field()
        plot_Density()
        plot_E(i)
        fig.savefig('figs/im' + str(int(i/figStride)) + '.png', dpi=300)
    sim.advance(timeStep=timeStep)