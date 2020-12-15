#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu as lu

meshManager = dg.MeshManager()
meshManager.readMesh('./input/box.msh')

# Numerical parameters:
N = 1

nodes = dg.TriangleNodesProvisioner(N, meshManager)
ctx = nodes.dgContext()
x = ctx.x
y = ctx.y

poisson2d = dg.Poisson2DSparseMatrix(ctx, meshManager)

MM = poisson2d.getMM()
OP = poisson2d.getOP()

# create scipy sparse mats
dof = ctx.numLocalPoints*ctx.numElements
MMsp = csc_matrix((MM[:, 2], (MM[:, 0], MM[:, 1])), shape=(dof, dof), dtype=np.float)

rows = OP[:, 0]
cols = OP[:, 1]
vals = OP[:, 2]

# bordering trick - ensures Pure Neuman problem is solvable.
for j in range(0, dof):
    rows = np.append(rows, [dof])
    cols = np.append(cols, [j])
    vals = np.append(vals, [1.0])

    cols = np.append(cols, [dof])
    rows = np.append(rows, [j])
    vals = np.append(vals, [1.0])

cols = np.append(cols, [dof])
rows = np.append(rows, [dof])
vals = np.append(vals, [0.0])

# Pass sparse triplet into csc_matrix constructor.
OPsp = csc_matrix((vals, (rows, cols)), shape=(dof + 1, dof + 1))

# factorize
linSolver = lu(-OPsp)

# Define right-hand side; compute load vector (mass matrix times forcing).
rhs1 = np.sin(0.5*np.pi*x, dtype=np.float , order='C') * np.sin(0.5*np.pi*y, dtype=np.float , order='C')
rhs = rhs1.flatten('F')
rhs = MMsp.dot(rhs)
rhs = np.append(rhs, [0.0])

# solve x = A \ b
soln = linSolver.solve(rhs)
soln = np.reshape(soln[:-1], (ctx.numLocalPoints, ctx.numElements), order='F')

# Output solution statistics.
print(f"max: {np.max(soln)}")
print(f"min: {np.min(soln)}")
print(f"mean: {np.sum(soln)/dof}")

fields = dict()
# Need to ensure fields send to the outputter are Row-Major ('C'-ordering),
# since blitzdg internally assumes they are.
fields["u"] = np.array(soln, dtype=np.float, order='C')
fields["f"] = np.array(rhs1, dtype=np.float, order='C')
outputter = dg.VtkOutputter(nodes)
outputter.writeFieldsToFiles(fields, 0)
