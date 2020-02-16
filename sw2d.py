#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
import matplotlib.pyplot as plt
from pprint import pprint

def sw2dComputeFluxes(h, hu, hv, hN, g, H):
    #h equation
    F1 = hu
    G1 = hv

    # Get velocity fields
    u = hu / h
    v = hv / h

    # hu equation

    F2 = hu*u + 0.5*g*h*h
    G2 = hu*v

    # hv equation
    F3 = G2
    G3 = hv*v + 0.5*g*h*h

    # N (tracer) equation
    F4 = hN*u
    G4 = hN*v

    return ((F1,F2,F3,F4),(G1,G2,G3,G4))

def sw2dComputeRHS(h, hu, hv, hN, g, H, f, ctx):
    vmapM = ctx.vmapM
    vmapP = ctx.vmapP
    BCmap = ctx.BCmap
    nx = ctx.nx
    ny = ctx.ny
    rx = ctx.rx
    sx = ctx.sx
    ry = ctx.ry
    sy = ctx.sy
    Dr = ctx.Dr
    Ds = ctx.Ds
    Nfp = ctx.numFacePoints

    Lift = ctx.Lift
    Fscale = ctx.Fscale

    hC = h.flatten('F')
    huC = hu.flatten('F')
    hvC = hv.flatten('F')
    hNC = hN.flatten('F')
    nxC = nx.flatten('F')
    nyC = ny.flatten('F')
    
    mapW = BCmap[3]

    # get field values along elemental faces.
    hM = hC[vmapM]
    hP = hC[vmapP]

    huM = huC[vmapM]
    huP = huC[vmapP]

    hvM = hvC[vmapM]
    hvP = hvC[vmapP]

    hNM = hNC[vmapM]
    hNP = hNC[vmapP]

    nxW = nxC[mapW]
    nyW = nyC[mapW]

    # set bc's (no normal flow thru the walls).
    huP[mapW] = huM[mapW] - 2*nxW*(huM[mapW]*nxW + hvM[mapW]*nyW)
    hvP[mapW] = hvM[mapW] - 2*nyW*(huM[mapW]*nxW + hvM[mapW]*nyW)

    # compute jump in states
    dh = hM - hP
    dhu = huM - huP
    dhv = hvM - hvP
    dhN = hNM - hNP

    ((F1M,F2M,F3M,F4M),(G1M,G2M,G3M,G4M)) = sw2dComputeFluxes(hM, huM, hvM, hNM, g, H)
    ((F1P,F2P,F3P,F4P),(G1P,G2P,G3P,G4P)) = sw2dComputeFluxes(hP, huP, hvP, hNP, g, H)
    ((F1,F2,F3,F4),(G1,G2,G3,G4)) = sw2dComputeFluxes(h, hu, hv, hN, g, H)

    uM = huM/hM 
    vM = hvM/hM

    uP = huP/hP
    vP = hvP/hP

    spdM = np.sqrt(uM*uM + vM*vM) + np.sqrt(g*hM)
    spdP = np.sqrt(uP*uP + vP*vP) + np.sqrt(g*hP)

    spdMax = np.max(np.array([spdM, spdP]), axis=0)

    # spdMax = np.max(spdMax)
    lam = np.reshape(spdMax, (ctx.numFacePoints, ctx.numFaces*ctx.numElements), order='F')
    lamMaxMat = np.outer(np.ones((Nfp, 1), dtype=np.float), np.max(lam, axis=0))
    spdMax = lamMaxMat.flatten('F')

    # strong form: Compute flux jump vector. (fluxM - numericalFlux ) dot n
    dFlux1 = 0.5*((F1M - F1P)*nxC + (G1M-G1P)*nyC - spdMax*dh)
    dFlux2 = 0.5*((F2M - F2P)*nxC + (G2M-G2P)*nyC - spdMax*dhu)
    dFlux3 = 0.5*((F3M - F3P)*nxC + (G3M-G3P)*nyC - spdMax*dhv)
    dFlux4 = 0.5*((F4M - F4P)*nxC + (G4M-G4P)*nyC - spdMax*dhN)

    dFlux1Mat = np.reshape(dFlux1, (Nfp*ctx.numFaces, K), order='F')
    dFlux2Mat = np.reshape(dFlux2, (Nfp*ctx.numFaces, K), order='F')
    dFlux3Mat = np.reshape(dFlux3, (Nfp*ctx.numFaces, K), order='F')
    dFlux4Mat = np.reshape(dFlux4, (Nfp*ctx.numFaces, K), order='F')

    # Flux divergence:
    RHS1 = -(rx*np.dot(Dr, F1) + sx*np.dot(Ds, F1))
    RHS1+= -(ry*np.dot(Dr, G1) + sy*np.dot(Ds, G1))

    RHS2 = -(rx*np.dot(Dr, F2) + sx*np.dot(Ds, F2))
    RHS2+= -(ry*np.dot(Dr, G2) + sy*np.dot(Ds, G2))

    RHS3 = -(rx*np.dot(Dr, F3) + sx*np.dot(Ds, F3))
    RHS3+= -(ry*np.dot(Dr, G3) + sy*np.dot(Ds, G3))

    RHS4 = -(rx*np.dot(Dr, F4) + sx*np.dot(Ds, F4))
    RHS4+= -(ry*np.dot(Dr, G4) + sy*np.dot(Ds, G4))

    surfaceRHS1 = Fscale*dFlux1Mat
    surfaceRHS2 = Fscale*dFlux2Mat
    surfaceRHS3 = Fscale*dFlux3Mat
    surfaceRHS4 = Fscale*dFlux4Mat

    RHS1 += np.dot(Lift, surfaceRHS1)
    RHS2 += np.dot(Lift, surfaceRHS2)
    RHS3 += np.dot(Lift, surfaceRHS3)
    RHS4 += np.dot(Lift, surfaceRHS4)

    # Add source terms
    RHS2 += f*hv
    RHS3 -= f*hu

    return (RHS1, RHS2, RHS3, RHS4)

# Main solver:
# set scaled density jump.
drho = 1.0025 - 1.000

# compute reduced gravity
g = drho*9.81

# set f-plane Coriolis frequency.
f = 7.88e-5

finalTime = 24.0*3600*3
numOuts = 200
t = 0.0

meshManager = dg.MeshManager()
meshManager.readMesh('input/R_8km_circle.msh')

# Numerical parameters:
NOrder = 6

filtOrder = 4
filtCutoff = 0.9*NOrder

nodes = dg.TriangleNodesProvisioner(NOrder, meshManager)
nodes.buildFilter(filtCutoff, filtOrder)

outputter = dg.VtkOutputter(nodes)

ctx = nodes.dgContext()

x = ctx.x
y = ctx.y

Np = ctx.numLocalPoints
K = ctx.numElements

Filt = ctx.filter


eta = -2.5*(x/8000.0)
u   = np.zeros([Np, K], dtype=np.float, order='C')
v   = np.zeros([Np, K], dtype=np.float, order='C')
H   = 10*np.ones([Np, K], dtype=np.float, order='C')
Nrad = 2e3
Nx = 2000.0
Ny = 2500.0
# N   = np.exp(-(((x-Nx)/Nrad)**2 + ((y-Ny)/Nrad)**2))
N   = np.exp(-((y-Ny)/Nrad)**2)

h = H + eta
hu = h*u
hv = h*v
hN = h*N

# setup fields dictionary for outputting.
fields = dict()
fields["eta"] = eta
fields["u"] = u
fields["v"] = v
fields["N"] = N
outputter.writeFieldsToFiles(fields, 0)

c = np.sqrt(g*h)
CFL = 0.8
dt = CFL / np.max( ((NOrder+1)**2)*0.5*np.abs(ctx.Fscale.flatten('F'))*(c.flatten('F')[ctx.vmapM]  + np.sqrt(((u.flatten('F'))[ctx.vmapM])**2 + ((v.flatten('F'))[ctx.vmapM])**2)))

numSteps = int(np.ceil(finalTime/dt))
outputInterval = int(numSteps / numOuts)

step = 0
while t < finalTime:

    (RHS1,RHS2,RHS3,RHS4) = sw2dComputeRHS(h, hu, hv, hN, g, H, f, ctx)

    RHS1 = np.dot(Filt, RHS1)
    RHS2 = np.dot(Filt, RHS2)
    RHS3 = np.dot(Filt, RHS3)
    RHS4 = np.dot(Filt, RHS4)
    
    # predictor
    h1  = h + 0.5*dt*RHS1
    hu1 = hu + 0.5*dt*RHS2
    hv1 = hv + 0.5*dt*RHS3
    hN1 = hN + 0.5*dt*RHS4

    (RHS1,RHS2,RHS3,RHS4) = sw2dComputeRHS(h1, hu1, hv1, hN1, g, H, f, ctx)

    RHS1 = np.dot(Filt, RHS1)
    RHS2 = np.dot(Filt, RHS2)
    RHS3 = np.dot(Filt, RHS3)
    RHS4 = np.dot(Filt, RHS4)

    # corrector - Update solution
    h += dt*RHS1
    hu += dt*RHS2
    hv += dt*RHS3
    hN += dt*RHS4

    h_max = np.max(np.abs(h))
    if h_max > 1e8  or np.isnan(h_max):
        raise Exception("A numerical instability has occurred.")

    t += dt
    step += 1

    if step % outputInterval == 0 or step == numSteps:
        print('Outputting at t=' + str(t))
        eta = h-H
        fields["eta"] = eta
        fields["u"] = hu/h
        fields["v"] = hv/h
        fields["N"] = hN/h
        outputter.writeFieldsToFiles(fields, step)
