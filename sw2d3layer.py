#!/usr/bin/python3
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg

def sw3lrComputeHyperbolicFluxes(h1, h1u1, h1v1, h2, h2u2, h2v2, h3, h3u3, h3v3, g):
    # h1 equation
    F1 = h1u1
    G1 = h1v1

    # h1u1 equation
    F2 = (h1u1*h1u1)/h1 + 0.5*g*h1*h1
    G2 = (h1u1*h1v1)/h1

    # h1v1 equation
    F3 = G2
    G3 = (h1v1*h1v1)/h1 + 0.5*g*h1*h1

    # h2 equation
    F4 = h2u2
    G4 = h2v2

    # h2u2 equation
    F5 = (h2u2*h2u2)/h2 + 0.5*g*h2*h2
    G5 = (h2u2*h2v2)/h2

    # h2v2 equation
    F6 = G5
    G6 = (h2v2*h2v2)/h2 + 0.5*g*h2*h2

    # h3 equation
    F7 = h3u3
    G7 = h3v3

    # h3u3 equation
    F8 = (h3u3*h3u3)/h3 + 0.5*g*h3*h3
    G8 = (h3u3*h3v3)/h3

    # h3v3 equation
    F9 = G8
    G9 = (h3v3*h3v3)/h3 + 0.5*g*h3*h3


    return ((F1,F2,F3,F4,F5,F6,F7,F8,F9)
        ,(G1,G2,G3,G4,G5,G6,G7,G8,G9))

def sw2dComputeRHS(h1, h1u1, h1v1, h2, h2u2, h2v2, h3, h3u3, h3v3, g, Bx, By, ctx, r):
    rx = ctx.rx
    sx = ctx.sx
    ry = ctx.ry
    sy = ctx.sy
    Dr = ctx.Dr
    Ds = ctx.Ds
    vmapM = ctx.vmapM

    mapW = ctx.BCmap[3]
    vmapW = vmapM[mapW]

    ((F1,F2,F3,F4,F5,F6,F7,F8,F9),
        (G1,G2,G3,G4,G5,G6,G7,G8,G9)) = sw3lrComputeHyperbolicFluxes(h1, h1u1, h1v1, 
            h2, h2u2, h2v2, h3, h3u3, h3v3, g)

    # Compute (Hyperbolic) Flux divergence, and add source terms if any.
    RHS1 = -(rx*np.dot(Dr, F1) + sx*np.dot(Ds, F1))
    RHS1+= -(ry*np.dot(Dr, G1) + sy*np.dot(Ds, G1))

    RHS2 = -(rx*np.dot(Dr, F2) + sx*np.dot(Ds, F2))
    RHS2+= -(ry*np.dot(Dr, G2) + sy*np.dot(Ds, G2))
    RHS2+= -g*h1*Bx

    RHS3 = -(rx*np.dot(Dr, F3) + sx*np.dot(Ds, F3))
    RHS3+= -(ry*np.dot(Dr, G3) + sy*np.dot(Ds, G3))
    RHS3+= -g*h1*By

    # layer 2:
    RHS4 = -(rx*np.dot(Dr, F4) + sx*np.dot(Ds, F4))
    RHS4+= -(ry*np.dot(Dr, G4) + sy*np.dot(Ds, G4))

    RHS5 = -(rx*np.dot(Dr, F5) + sx*np.dot(Ds, F5))
    RHS5+= -(ry*np.dot(Dr, G5) + sy*np.dot(Ds, G5))
    RHS5+= -g*h2*Bx - g*(2*r/(1+r))*h2*(rx*np.dot(Dr, h1) + sx*np.dot(Ds, h1))
    RHS5+= -g*h2*(rx*np.dot(Dr, h3) + sx*np.dot(Ds, h3))

    RHS6 = -(rx*np.dot(Dr, F6) + sx*np.dot(Ds, F6))
    RHS6+= -(ry*np.dot(Dr, G6) + sy*np.dot(Ds, G6))
    RHS6+= -g*h2*By - g*(2*r/(1+r))*h2*(ry*np.dot(Dr, h1) + sy*np.dot(Ds, h1))
    RHS6+= -g*h2*(ry*np.dot(Dr, h3) + sy*np.dot(Ds, h3))

    # layer 3:
    RHS7 = -(rx*np.dot(Dr, F7) + sx*np.dot(Ds, F7))
    RHS7+= -(ry*np.dot(Dr, G7) + sy*np.dot(Ds, G7))

    RHS8 = -(rx*np.dot(Dr, F8) + sx*np.dot(Ds, F8))
    RHS8+= -(ry*np.dot(Dr, G8) + sy*np.dot(Ds, G8))
    RHS8+= -g*h3*Bx - g*0.5*(1+r)*h3*(rx*np.dot(Dr, h3) + sx*np.dot(Ds, h3))
    RHS8+= -g*r*h3*(rx*np.dot(Dr, h1) + sx*np.dot(Ds, h1))

    RHS9 = -(rx*np.dot(Dr, F9) + sx*np.dot(Ds, F9))
    RHS9+= -(ry*np.dot(Dr, G9) + sy*np.dot(Ds, G9))
    RHS9+= -g*h3*By - g*0.5*(1+r)*h3*(ry*np.dot(Dr, h3) + sy*np.dot(Ds, h3))
    RHS9+= -g*r*h3*(ry*np.dot(Dr, h1) + sy*np.dot(Ds, h1))

    return (RHS1, RHS2, RHS3)

# Main solver:
if __name__ == '__main__':
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
