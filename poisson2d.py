#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
import matplotlib.pyplot as plt
from pprint import pprint
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
MMsp = csc_matrix((MM[:, 2], (MM[:, 0], MM[:, 1])), shape=(dof, dof))
OPsp = csc_matrix((OP[:, 2], (OP[:, 0], OP[:, 1])), shape=(dof, dof))


# factorize
linSolver = lu(OPsp)

# solve x = A \ b

rhs = np.cos(x)*np.cos(y)
rhs = rhs.flatten('F')

soln = linSolver.solve(rhs)

print(soln)

plt.figure()
plt.spy(OPsp)
plt.show()

outputter = dg.VtkOutputter(nodes)


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