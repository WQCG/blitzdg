#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from lib import pyblitzdg as dg
import matplotlib.pyplot as plt
from pprint import pprint

def sw2dComputeFluxes(h, hu, hv, g, H):
    #h equation
    F1 = hu
    G1 = hv

    # hu equation
    F2 = (hu*hu)/h + 0.5*g*h*h
    G2 = (hu*hv)/h

    # hv equation
    F3 = G2
    G3 = (hv*hv)/h + 0.5*g*h*h

    return ((F1,F2,F3),(G1,G2,G3))

def sw2dComputeRHS(h, hu, hv, g, H, ctx):
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

    nxW = nxC[mapW]
    nyW = nyC[mapW]

    # set bc's (no normal flow thru the walls).
    huP[mapW] = huM[mapW] - 2*nxW*(huM[mapW]*nxW + hvM[mapW]*nyW)
    hvP[mapW] = hvM[mapW] - 2*nyW*(huM[mapW]*nxW + hvM[mapW]*nyW)

    # compute jump in states
    dh = hM - hP
    dhu = huM - huP
    dhv = hvM - hvP

    ((F1M,F2M,F3M),(G1M,G2M,G3M)) = sw2dComputeFluxes(hM, huM, hvM, g, H)
    ((F1P,F2P,F3P),(G1P,G2P,G3P)) = sw2dComputeFluxes(hP, huP, hvP, g, H)
    ((F1,F2,F3),(G1,G2,G3)) = sw2dComputeFluxes(h, hu, hv, g, H)

    uM = huM/hM 
    vM = hvM/hM

    uP = huP/hP
    vP = hvP/hP

    spdM = np.sqrt(uM*uM + vM*vM) + np.sqrt(g*hM)
    spdP = np.sqrt(uP*uP + vP*vP) + np.sqrt(g*hP)

    spdMax = np.max(np.array([spdM, spdP]), axis=0)

    # spdMax = np.max(spdMax)
    lam = np.reshape(spdMax, (ctx.numFacePoints, ctx.numFaces*ctx.numElements), order='F')
    lamMaxMat = np.outer(np.ones((Nfp, 1), dtype=np.dtype('Float64')), np.max(lam, axis=0))
    spdMax = lamMaxMat.flatten('F')

    # strong form: Compute flux jump vector. (fluxM - numericalFlux ) dot n
    dFlux1 = 0.5*((F1M - F1P)*nxC + (G1M-G1P)*nyC - spdMax*dh)
    dFlux2 = 0.5*((F2M - F2P)*nxC + (G2M-G2P)*nyC - spdMax*dhu)
    dFlux3 = 0.5*((F3M - F3P)*nxC + (G3M-G3P)*nyC - spdMax*dhv)

    dFlux1Mat = np.reshape(dFlux1, (Nfp*ctx.numFaces, K), order='F')
    dFlux2Mat = np.reshape(dFlux2, (Nfp*ctx.numFaces, K), order='F')
    dFlux3Mat = np.reshape(dFlux3, (Nfp*ctx.numFaces, K), order='F')

    # Flux divergence:
    RHS1 = -(rx*np.dot(Dr, F1) + sx*np.dot(Ds, F1))
    RHS1+= -(ry*np.dot(Dr, G1) + sy*np.dot(Ds, G1))

    RHS2 = -(rx*np.dot(Dr, F2) + sx*np.dot(Ds, F2))
    RHS2+= -(ry*np.dot(Dr, G2) + sy*np.dot(Ds, G2))

    RHS3 = -(rx*np.dot(Dr, F3) + sx*np.dot(Ds, F3))
    RHS3+= -(ry*np.dot(Dr, G3) + sy*np.dot(Ds, G3))

    # to check are rx,ry,sx,sy the same as cpp? what about fscale?
    # then check dflux mat's and compare. I think lift is the same but can check again.

    surfaceRHS1 = Fscale*dFlux1Mat
    surfaceRHS2 = Fscale*dFlux2Mat
    surfaceRHS3 = Fscale*dFlux3Mat

    RHS1 += np.dot(Lift, surfaceRHS1)
    RHS2 += np.dot(Lift, surfaceRHS2)
    RHS3 += np.dot(Lift, surfaceRHS3)

    return (RHS1, RHS2, RHS3)

# Main solver:
g = 9.81
finalTime = 40.0
t = 0.0

meshManager = dg.MeshManager()
meshManager.readMesh('./input/box.msh')

# Numerical parameters:
N = 1
CFL = 0.45

filtOrder = 4
filtCutoff = 0.9*N

nodes = dg.TriangleNodesProvisioner(N, meshManager)
nodes.buildFilter(filtCutoff, filtOrder)

outputter = dg.VtkOutputter(nodes)

ctx = nodes.dgContext()

x = ctx.x
y = ctx.y


Np = ctx.numLocalPoints
K = ctx.numElements

Filt = ctx.filter

eta = np.exp(-10*(x*x) -10*(y*y), dtype=np.dtype('Float64') , order='C')
#eta = -1*(x/1500.0)
u   = np.zeros([Np, K], dtype=np.dtype('Float64'), order='C')
v   = np.zeros([Np, K], dtype=np.dtype('Float64'), order='C')
H   = 10*np.ones([Np, K], dtype=np.dtype('Float64'), order='C')

h = H + eta
hu = h*u
hv = h*v

# setup fields dictionary for outputting.
fields = dict()
fields["eta"] = eta
fields["u"] = u
fields["v"] = v
outputter.writeFieldsToFiles(fields, 0)

c = np.sqrt(g*h)
dt =0.000724295
#dt = CFL*dx/np.max(abs(c))

step = 0
while t < finalTime:

    
    (RHS1,RHS2,RHS3) = sw2dComputeRHS(h, hu, hv, g, H, ctx)
    
    # predictor
    h1  = h + 0.5*dt*RHS1
    hu1 = hu + 0.5*dt*RHS2
    hv1 = hv + 0.5*dt*RHS3

    (RHS1,RHS2,RHS3) = sw2dComputeRHS(h1, hu1, hv1, g, H, ctx)

    # corrector - Update solution
    h += dt*RHS1
    hu += dt*RHS2
    hv += dt*RHS3

    h_max = np.max(np.abs(h))
    if h_max > 1e8  or np.isnan(h_max):
        raise Exception("A numerical instability has occurred.")

    t += dt
    step += 1

    print('t=' + str(t))

    eta = h-H
    fields["eta"] = eta
    fields["u"] = hu/h
    fields["v"] = hv/h
    outputter.writeFieldsToFiles(fields, step)
