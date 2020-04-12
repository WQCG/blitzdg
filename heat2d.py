#!/usr/bin/python3
'''
Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import splu as lu
import cghelpers as cg

meshManager = dg.MeshManager()
meshManager.readMesh('./input/trapezoid.msh')

Verts = meshManager.vertices
EToV = meshManager.elements
bcType = meshManager.bcType


# Correct boundary condition table for the different types
# TODO, store magic integers like 6 - Dirichlet and 7 - Neuman in an enum prop.
faces = np.array([[0, 1], [1, 2], [2, 0]])
for e in range(0, len(bcType)):
    for (faceInd, f) in enumerate(faces):
        v1 = EToV[e, f[0]]
        v2 = EToV[e, f[1]]
        v1x = Verts[v1, 0]
        v1y = Verts[v1, 1]
        v2x = Verts[v2, 0]
        v2y = Verts[v2, 1]

        midx = 0.5*(v1x + v2x)
        midy = 0.5*(v1y + v2y)

        if  np.abs(midy - 1.3) < 1e-12:
            bcType[e, faceInd] = 6
        elif np.abs(midy) < 1e-12:
            bcType[e, faceInd] = 7   # 7 - bottom neuman
        elif np.abs(midy - 1.3*midx) < 1e-12:
            bcType[e, faceInd] = 6
        elif np.abs(midy + 1.3*midx - 26/5) < 1e-12:
            bcType[e, faceInd] = 6

# Update BC table in the backend
meshManager.setBCType(bcType)

#set diffusivity param
kappa = 1
#kappa = 1.0 - 0.9*np.exp(-100*(x-2.0)**2)

# Numerical parameters:
N = 3

nodes = dg.TriangleNodesProvisioner(N, meshManager)
ctx = nodes.dgContext()
x = ctx.x
y = ctx.y
outputter = dg.VtkOutputter(nodes)

CFL = 0.85

dt = CFL / np.max( ((N+1)**2)*0.5*(np.abs(ctx.Fscale.flatten('F'))**2)*kappa)

Np = ctx.numLocalPoints
K = ctx.numElements
Nfaces = ctx.numFaces
Nfp = ctx.numFacePoints

# identify boundary face nodes
vmapM = ctx.vmapM

bcMap = ctx.BCmap
mapD = bcMap[6]
mapN = bcMap[7]

Fx = x.flatten('F')[vmapM]
Fy = y.flatten('F')[vmapM]

bottomLeft = np.where(np.abs(Fx) < 1e-12)[0]

top  = np.where(np.abs(Fy - 1.3) < 1e-12)[0]
left = np.where(np.abs(Fy - 1.3*Fx) < 1e-12)[0]
right = np.where(np.abs(Fy + 1.3*Fx - 26/5) < 1e-12)[0]

bottom = np.where(np.abs(Fy) < 1e-12)[0] #chop first and last points from neuman bc
bottomNoRight = np.where(Fx[bottom] < 4)
bottom = bottom[bottomNoRight]

bcType = meshManager.bcType

EToV = meshManager.elements
Verts = meshManager.vertices

# Constructing poisson2D sparse matrix
poisson2d = dg.Poisson2DSparseMatrix(ctx, meshManager)

#Neuman - homogeneous
qbc = np.zeros((Nfp*Nfaces*K))

# Dirichlet BCs
ubc = np.zeros((Nfp*Nfaces*K))
ubc[right] = np.ones(right.shape) - (10/13)*Fy[right]
ubc[top]   = 0.25*(np.exp(4*(Fx[top]-3)) - 1)
ubc[left]  = -0.25*np.ones(left.shape)

# Put in the shape required for the backend.
qbc = np.reshape(qbc, (Nfp*Nfaces, K), order='F')
ubc = np.reshape(ubc, (Nfp*Nfaces, K), order='F')

# Compute RHS
bcRhs = poisson2d.buildBcRhs(ctx, meshManager, np.array(ubc, dtype=np.float64, order='C'), np.array(qbc, dtype=np.float64, order='C'))

# get sparse triplet global operators
MM = poisson2d.getMM()
OP = poisson2d.getOP()

# create scipy sparse matrices.
dof = ctx.numLocalPoints*ctx.numElements
MMsp = csc_matrix((MM[:, 2], (MM[:, 0], MM[:, 1])), shape=(dof, dof), dtype=np.float64)

rows = OP[:, 0]
cols = OP[:, 1]
vals = OP[:, 2]

# scale laplacian assuming const dt.
vals = kappa*dt*(vals)

# Pass sparse triplet into csc_matrix constructor.
scaledLap = csc_matrix((vals, (rows, cols)), shape=(dof, dof))

# form full LHS operator.
OP = MMsp + scaledLap

# factorize
linSolver = lu(OP)

# Define initial condition.
p = np.zeros((Np, K))

bcContrib = bcRhs.flatten('F')

fields = dict()
fields["p"] = np.array(p, dtype=np.float64, order='C')
outputter.writeFieldsToFiles(fields, 0)

for j in range(1, 2000):

    # compute load vector (MassMatrix times p)
    rhs = p.flatten('F')
    MMrhs = MMsp.dot(rhs)

    # solve x = A \ b
    soln = linSolver.solve(MMrhs + bcContrib)
    soln = np.reshape(soln, (ctx.numLocalPoints, ctx.numElements), order='F')

    # output every 10th timestep
    if j % 10 == 0:
        print("Writing output at t=" + str(dt*j))

        # Need to ensure fields send to the outputter are Row-Major ('C'-ordering),
        # since blitzdg internally assumes they are.
        fields["p"] = np.array(soln/np.max(np.abs(soln)), dtype=np.float64, order='C')
        outputter.writeFieldsToFiles(fields, j)

    # rotate soln and initial state for next timestep
    p = soln