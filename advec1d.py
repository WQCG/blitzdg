'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from lib import pyblitzdg
import matplotlib.pyplot as plt

def advec1dComputeRHS(u,c, nodes1d):
    vmapM = nodes1d.vmapM
    vmapP = nodes1d.vmapP
    mapI = nodes1d.mapI
    mapO = nodes1d.mapO
    nx = nodes1d.nx
    rx = nodes1d.rx
    Dr = nodes1d.Dr
    Lift = nodes1d.Lift
    Fscale = nodes1d.Fscale

    uC = u.flatten('F')
    nxC = nx.flatten('F')
    RHS = 0*vmapM
    uM = uC[vmapM]
    uP = uC[vmapP]

    # set bc's (outflow condition, and no inflow).
    uP[mapO] = uM[mapO]
    uP[mapI] = 0

    # compute jump in flux with Lax-Friedrichs numerical flux.
    alpha = 0
    du = (uM - uP)*0.5*(c*nxC - (1-alpha)*np.abs(c*nxC)); 

    RHS =-c*rx*np.dot(Dr, u)
    RHS += np.dot(Lift, Fscale*np.reshape(du, nx.shape, 'F'))

    return RHS

# Main solver:

xmin =-1.0
xmax = 4.0
c = 0.1

finalTime = 20.0
t = 0.0

# Numerical parameters:
N = 1
K = 60
CFL = 0.25

nodes1d = pyblitzdg.Nodes1DProvisioner(N, K, xmin, xmax)
nodes1d.buildNodes()
nodes1d.computeJacobian()
x = nodes1d.xGrid

u = np.exp(-10*(x*x))
xplot = np.squeeze(np.reshape(x, ((N+1)*K, 1)))
uplot = np.squeeze(np.reshape(u, ((N+1)*K, 1)))

xplot = np.squeeze(np.reshape(x, ((N+1)*K, 1)))
uplot = np.squeeze(np.reshape(u, ((N+1)*K, 1)))

min_dx = x[1,0] - x[0,0]

dt = CFL*min_dx/abs(c)

plt.figure()
plt.ion()
plt.show()
while t < finalTime:

    # Calculate Right-hand side.
    RHS = advec1dComputeRHS(u,c, nodes1d)

    # Update solution.
    u += dt*RHS

    u_max = np.max(np.max(np.abs(u)))
    if u_max > 1e8  or np.isnan(u_max):
        raise("A numerical instability has occurred!")

    t += dt

    xplot = x.flatten('F')
    uplot = u.flatten('F')

    plt.clf()
    plt.plot(x, u, '.-b')
    plt.draw()
    plt.pause(1.e-4)

