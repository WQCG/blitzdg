#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix, eye
from scipy.sparse.linalg import splu as lu

meshManager = dg.MeshManager()
meshManager.readMesh('./input/box.msh')

# Numerical parameters:
N = 1

nodes = dg.TriangleNodesProvisioner(N, meshManager)
ctx = nodes.dgContext()
x = ctx.x
y = ctx.y

kappa = 0.05
dt = 0.005

poisson2d = dg.Poisson2DSparseMatrix(ctx, meshManager)

MM = poisson2d.getMM()
OP = poisson2d.getOP()

# create scipy sparse mats
dof = ctx.numLocalPoints*ctx.numElements
MMsp = csc_matrix((MM[:, 2], (MM[:, 0], MM[:, 1])), shape=(dof, dof), dtype=np.dtype('Float64'))

rows = OP[:, 0]
cols = OP[:, 1]
vals = OP[:, 2]

vals = kappa*dt*(-vals)

# Pass sparse triplet into csc_matrix constructor.
scaledLap = csc_matrix((vals, (rows, cols)), shape=(dof, dof))

OPsp = MMsp - scaledLap

# factorize
linSolver = lu(OPsp)

# Define initial condition.
u = 2.5*np.exp(-((x/0.1)**2.0 + (y/0.1)**2.0), dtype=np.dtype('Float64'), order='C')

for j in range(0, 500):

    # compute load vector (MassMatrix times u)
    rhs = u.flatten('F')
    rhs = MMsp.dot(rhs)

    # solve x = A \ b
    soln = linSolver.solve(rhs)
    soln = np.reshape(soln, (ctx.numLocalPoints, ctx.numElements), order='F')

    # Output solution statistics.
    #print(f"max: {np.max(soln)}")
    #print(f"min: {np.min(soln)}")
    #print(f"mean: {np.sum(soln)/dof}")

    # output every 10th timestep
    if j % 10 == 0:
        fields = dict()
        # Need to ensure fields send to the outputter are Row-Major ('C'-ordering),
        # since blitzdg internally assumes they are.
        fields["u"] = np.array(soln, dtype=np.dtype('Float64'), order='C')
        outputter = dg.VtkOutputter(nodes)
        outputter.writeFieldsToFiles(fields, j)

    # rotate soln and initial state for next timestep
    u = soln