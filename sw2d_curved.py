#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
import pyblitzdg as dg
import matplotlib.pyplot as plt
from pprint import pprint
from scipy.interpolate import splev, splrep, interp1d, griddata

from swhelpers.flux import sw2dComputeFluxes
from swhelpers.maps import makeMapsPeriodic, correctBCTable
from swhelpers.rhs import sw2dComputeRHS_curved
from swhelpers.limiters import surfaceReconstruction, positivityPreservingLimiter2D

# Main solver:
# Set scaled density jump.
drho = 1.00100 - 1.000

# compute reduced gravity
g = drho*9.81

# Set depth.
H0 = 7.5
c0 = np.sqrt(g*H0)

finalTime = 24*3600
numOuts = 200
t = 0.0

meshManager = dg.MeshManager()
meshManager.readMesh('input/headlands_highres.msh')

Verts = meshManager.vertices
EToV = meshManager.elements
bcType = meshManager.bcType

# 2 = outflow.
bcType = correctBCTable(bcType, EToV, Verts, 2)
meshManager.setBCType(bcType)

# Numerical parameters:
NOrder = 4

filtOrder = 4
filtCutoff = 0.9*NOrder

nodes = dg.TriangleNodesProvisioner(NOrder, meshManager)
nodes.buildFilter(filtCutoff, filtOrder)
outputter = dg.VtkOutputter(nodes)
ctx = nodes.dgContext()

x = ctx.x
y = ctx.y

indN, indK = np.where(np.hypot(x, y) < 1.9e1)
centreIndN = indN[0]
centreIndK = indK[0]

xFlat = x.flatten('F')
yFlat = y.flatten('F')

mapW = ctx.BCmap[3]
vmapW = ctx.vmapM[mapW]
xW = xFlat[vmapW]
yW = yFlat[vmapW]

topInds = np.logical_and(yW > 200, np.logical_and(xW > 3250, xW < 4750))

xtop = xW[topInds]
ytop = yW[topInds]

isort = np.argsort(xtop)
xtop2 = xtop[isort]
ytop2 = ytop[isort]

# build arc-length parameter for headlands curve. 
s = [0.0]
for i in range(1, len(xtop2)):
    d = np.hypot(xtop2[i]-xtop2[i-1], ytop2[i]-ytop2[i-1])
    s.append(s[i-1] + d)

s = np.array(s)

s256 = np.linspace(s[0], s[-1], 128)
ss = np.linspace(s[0], s[-1], 4096)

xx = griddata(s, xtop2, s256, 'linear')
yy = griddata(s, ytop2, s256, 'linear')

# build parametric curve of 'top curve'/headlands curve.
splx = splrep(s256, xx)
sply = splrep(s256, yy)

xTopSmooth = splev(ss, splx, ext=2)
yTopSmooth = splev(ss, sply, ext=2)

bcInds = np.where(bcType.flatten('F') > 0)
bcFaces = np.transpose(np.unravel_index(bcInds, (ctx.numElements, ctx.numFaces), order='F'))

modifiedVerts = np.zeros(Verts.shape[0])

curvedFaces = []
for bcFace in bcFaces:

    el = bcFace[0][0]
    face = bcFace[0][1]
    # Get vertices along face
    v1ind = EToV[el, face]
    v2ind = EToV[el, (face+1) % ctx.numFaces]

    v1 = Verts[v1ind, :]
    v2 = Verts[v2ind, :]

    hyps1 = np.hypot(xTopSmooth - v1[0], yTopSmooth-v1[1])
    hyps2 = np.hypot(xTopSmooth - v2[0], yTopSmooth-v2[1])

    minInd1 = np.argmin(hyps1)
    minInd2 = np.argmin(hyps2)

    mytol = 200

    if hyps1[minInd1] > mytol or hyps2[minInd2] > mytol:
        # nothing to do here
        continue

    curvedFaces.append([el, face])

    # set new vertex coordinates to closest spline point coordinates
    newx1 = xTopSmooth[minInd1] 
    newy1 = yTopSmooth[minInd1]

    newx2 = xTopSmooth[minInd2]
    newy2 = yTopSmooth[minInd2]

    # update mesh vertex locations
    Verts[v1ind, 0] = newx1
    Verts[v1ind, 1] = newy1

    Verts[v2ind, 0] = newx2
    Verts[v2ind, 1] = newy2

    modifiedVerts[v1ind] = 1  
    modifiedVerts[v2ind] = 1


curvedEls = []
for face in curvedFaces:
    k = face[0]
    f = face[1]

    curvedEls.append(k)

    if f==0:
       v1 = EToV[k, 0]
       v2 = EToV[k, 1]
       vr = ctx.r
    elif f==1:
        v1 = EToV[k, 1]
        v2 = EToV[k, 2]
        vr = ctx.s
    elif f==2:
        v1 = EToV[k, 0]
        v2 = EToV[k, 2]
        vr = ctx.s

    fmsk = ctx.Fmask[:, f]

    fr = vr[ctx.Fmask[:, f]]
    x1 = Verts[v1, 0]
    y1 = Verts[v1, 1]
    x2 = Verts[v2, 0]
    y2 = Verts[v2, 1]

    v1_dists2 = (x1-xTopSmooth)**2 + (y1-yTopSmooth)**2
    v2_dists2 = (x2-xTopSmooth)**2 + (y2-yTopSmooth)**2

    v1s_inds = np.where(np.sqrt(v1_dists2) < 1.0e-8)
    v2s_inds = np.where(np.sqrt(v2_dists2) < 1.0e-8)

    if len(v1s_inds) == 0 or len(v2s_inds) == 0:
        continue

    if len(v1s_inds) > 1 and len(v2s_inds) > 1:
        raise Exception('bad parameterization.')

    # found degenerate interval.
    if v1s_inds[0] == v2s_inds[0]:
        continue
    
    if len(v1s_inds) == 1 and len(v2s_inds) == 1:
        t1 = ss[v1s_inds[0]] # set end-points of parameter-space interval.
        t2 = ss[v2s_inds[0]]
    else:
        raise Exception("Error tabulating spline interval end-points.")

    tLGL = 0.5*t1*(1-fr) + 0.5*t2*(1+fr)

    # Basically xnew - xold
    # where xnew is evaluated using the parameterization
    splxev = splev(tLGL, splx, ext=3)
    splyev = splev(tLGL, sply, ext=3) 
    x1 = x[fmsk, k]
    y1 = y[fmsk, k]
    fdx = splxev - x1
    fdy = splyev - y1

    # build 1D Vandermonde matrix for face nodes and volume nodes
    vand = dg.VandermondeBuilder()
    Vface, Vfinv = vand.buildVandermondeMatrix(fr, True, NOrder)
    Vvol,  = vand.buildVandermondeMatrix(vr, False, NOrder)

    # compute unblended volume deformations
    vdx = np.dot(Vvol, np.dot(Vfinv, fdx))
    vdy = np.dot(Vvol, np.dot(Vfinv, fdy))

    # compute blending functions for Gordon-Hall blending.
    r = ctx.r
    s = ctx.s
    ids = np.where(np.abs(1-vr) > 1.e-7)[0]
    blend = np.zeros(ids.shape)
    if f==0: blend = -(r[ids]+s[ids])/(1-vr[ids])
    if f==1: blend = (r[ids]+1)/(1-vr[ids])
    if f==2: blend = -(r[ids]+s[ids])/(1-vr[ids])
    
    # blend deformation to volume interior
    x[ids, k] += blend*vdx[ids]
    y[ids, k] += blend*vdy[ids]

print("Updating physical coordinates")
nodes.setCoordinates(x, y)

xr = np.dot(ctx.Dr, x)
yr = np.dot(ctx.Dr, y)
xs = np.dot(ctx.Ds, x)
ys = np.dot(ctx.Ds, y)

J = xr*ys - xs*yr
gauss_ctx = nodes.buildGaussFaceNodes(2*(NOrder+1))

BCmap = ctx.BCmap
mapW = ctx.BCmap[3]
mapO = ctx.BCmap[2]

vmapO = ctx.vmapM[mapO]

vmapW = ctx.vmapM[mapW]

xFlat = x.flatten('F')
yFlat = y.flatten('F')

gxFlat = np.dot(gauss_ctx.Interp, x).flatten('F')
gyFlat = np.dot(gauss_ctx.Interp, y).flatten('F')

gmapO = np.array(gauss_ctx.BCmap[2])

xO = xFlat[vmapO]
yO = yFlat[vmapO]

gxO = gxFlat[gmapO]
gyO = gyFlat[gmapO]

vmapM = ctx.vmapM
vmapP = ctx.vmapP

gmapM = gauss_ctx.mapM
gmapP = gauss_ctx.mapP

vmapM, vmapP = makeMapsPeriodic(vmapM, vmapP, vmapO, xFlat, yFlat, xO, yO)
gmapM, gmapP = makeMapsPeriodic(gmapM, gmapP, gmapO, gxFlat, gyFlat, gxO, gyO)

print("building cubature context")
cub_ctx = nodes.buildCubatureVolumeMesh(3*(NOrder+1))
print("done.")

Np = ctx.numLocalPoints
K = ctx.numElements
Filt = ctx.filter

f=7.8825e-5
amp = 0.5*.065*H0  # had 0.065*H0
L = 500
W = 200
eta = amp*np.exp(-((y-L)/W)**2)
etay = -2*amp*(y-L)*np.exp(-((y-L)/W)**2) / W**2
u = (-g/f)*etay

vort = -(ctx.ry*np.dot(ctx.Dr, u) + ctx.sy*np.dot(ctx.Ds, u))

# detect damping regions, and use to determine drag coefficient CD

xW = xFlat[vmapW]
yW = yFlat[vmapW]

CD = np.zeros((Np, K))

CD_max = 2.5e-3
CDflat = np.zeros((Np, K)).flatten('F')
length_tol = 2.5e2

for i, _ in enumerate(xFlat):
    xi = xFlat[i]
    yi = yFlat[i]

    dists = np.hypot(xi - xW, yi - yW)
    min_dist = np.min(dists)

    CDflat[i] = CD_max*0.5*(1-np.tanh((min_dist - 0.5*length_tol)/(0.1*length_tol)))

CD = np.array(np.reshape(CDflat, (Np, K), order='F'), order='C')

umax = np.max(u.flatten('F')) 
print("umax: ", umax)
print("rad:" , c0/f)
print("froude: ", umax/c0)
v = 0*eta

r = np.sqrt(x*x + y*y)

H = H0*np.ones((Np,K))
Dr = ctx.Dr 
Ds = ctx.Ds
rx = ctx.rx
ry = ctx.ry
sx = ctx.sx 
sy = ctx.sy

z = -H
zx = (rx*np.dot(Dr, z) + sx*np.dot(Ds, z))
zy = (ry*np.dot(Dr, z) + sy*np.dot(Ds, z))

Nrad = 3e2
Nx = 4000.0
Ny = 350.0
N   = np.exp(-(((x-Nx)/Nrad)**2 + ((y-Ny)/Nrad)**2))

h = H + eta
#h = 5*(0.5*(1 - np.tanh(x/Nrad)))
h[h < 1e-3] = 1e-3
#H = h - eta
eta = h - H
hu = h*u
hv = h*v
hN = h*N

# setup fields dictionary for outputting.
fields = dict()
fields["eta"] = eta
fields["u"] = u
fields["v"] = v
fields["B"] = N
fields["h"] = h
fields["vort"] = vort
outputter.writeFieldsToFiles(fields, 0)

Hbar = np.mean(H)
c = np.sqrt(g*Hbar)*np.ones((Np, K))
CFL = 0.75
dt = CFL / np.max( ((NOrder+1)**2)*0.5*np.abs(ctx.Fscale.flatten('F'))*(c.flatten('F')[vmapM]  + np.sqrt(((u.flatten('F'))[vmapM])**2 + ((v.flatten('F'))[vmapM])**2)))
print("dt=", str(dt))
# dt = 1.1
#dt = 0.05

numSteps = int(np.ceil(finalTime/dt))
outputInterval = 50


step = 0
print("Entering main time-loop")
fields = dict()
while t < finalTime:

    (RHS1,RHS2,RHS3,RHS4) = sw2dComputeRHS_curved(h, hu, hv, hN, zx, zy, g, H, f, CD, ctx, cub_ctx, gauss_ctx, curvedEls, J, gmapM, gmapP)
    
    RHS1 = np.dot(Filt, RHS1)
    RHS2 = np.dot(Filt, RHS2)
    RHS3 = np.dot(Filt, RHS3)
    RHS4 = np.dot(Filt, RHS4)
    
    # predictor
    h1  = h + 0.5*dt*RHS1
    hu1 = hu + 0.5*dt*RHS2
    hv1 = hv + 0.5*dt*RHS3
    hN1 = hN + 0.5*dt*RHS4

    #h1, hu1, hv1 = positivityPreservingLimiter2D(h1, hu1, hv1)
    #h1[h1 < 1e-3] = 1e-3

    (RHS1,RHS2,RHS3,RHS4) = sw2dComputeRHS_curved(h1, hu1, hv1, hN1, zx, zy, g, H, f, CD, ctx, cub_ctx, gauss_ctx, curvedEls, J, gmapM, gmapP)

    RHS1 = np.dot(Filt, RHS1)
    RHS2 = np.dot(Filt, RHS2)
    RHS3 = np.dot(Filt, RHS3)
    RHS4 = np.dot(Filt, RHS4)

    # corrector - Update solution
    h += dt*RHS1
    hu += dt*RHS2
    hv += dt*RHS3
    hN += dt*RHS4

    t += dt

    u = hu / h
    v = hv / h
    dt = CFL / np.max( ((NOrder+1)**2)*0.5*np.abs(ctx.Fscale.flatten('F'))*(c.flatten('F')[vmapM]  + np.sqrt(((u.flatten('F'))[vmapM])**2 + ((v.flatten('F'))[vmapM])**2)))


    h_max = np.max(np.abs(h))
    if h_max > 1e8  or np.isnan(h_max):
        raise Exception("A numerical instability has occurred.")

    step += 1

    if step % outputInterval == 0 or step == numSteps:
        print('Outputting at t=' + str(t))
        eta = h-H

        vort = ((ctx.rx*np.dot(ctx.Dr, v) + ctx.sx*np.dot(ctx.Ds, v)) - 
            (ctx.ry*np.dot(ctx.Dr, u) + ctx.sy*np.dot(ctx.Ds, u)))

        fields["eta"] = eta
        fields["u"] = hu/h
        fields["v"] = hv/h
        fields["B"] = hN/h
        fields["vort"] = vort
        outputter.writeFieldsToFiles(fields, step)
    
