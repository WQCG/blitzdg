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

def sw3lrComputeRHS(h1, h1u1, h1v1, h2, h2u2, h2v2, h3, h3u3, h3v3, g, Bx, By, ctx, r):
    rx = ctx.rx
    sx = ctx.sx
    ry = ctx.ry
    sy = ctx.sy
    Dr = ctx.Dr
    Ds = ctx.Ds
    vmapM = ctx.vmapM
    vmapP = ctx.vmapP
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
    RHS8+= -g*h3*Bx - g*0.5*(1+r)*h3*(rx*np.dot(Dr, h2) + sx*np.dot(Ds, h2))
    RHS8+= -g*r*h3*(rx*np.dot(Dr, h1) + sx*np.dot(Ds, h1))

    RHS9 = -(rx*np.dot(Dr, F9) + sx*np.dot(Ds, F9))
    RHS9+= -(ry*np.dot(Dr, G9) + sy*np.dot(Ds, G9))
    RHS9+= -g*h3*By - g*0.5*(1+r)*h3*(ry*np.dot(Dr, h2) + sy*np.dot(Ds, h2))
    RHS9+= -g*r*h3*(ry*np.dot(Dr, h1) + sy*np.dot(Ds, h1))

    # Bottom drag terms
    CD = 2.5e-3
    u3 = h3u3 / h3
    v3 = h3v3 / h3
    spd3 = np.sqrt(u3**2 + v3**2)

    RHS8 += -CD*h3u3*spd3
    RHS9 += -CD*h3v3*spd3

    # apply filter.
    Filt = ctx.filter
    RHS1 = np.dot(Filt, RHS1)
    RHS2 = np.dot(Filt, RHS2)
    RHS3 = np.dot(Filt, RHS3)
    RHS4 = np.dot(Filt, RHS4)
    RHS5 = np.dot(Filt, RHS5)
    RHS6 = np.dot(Filt, RHS6)
    RHS7 = np.dot(Filt, RHS7)
    RHS8 = np.dot(Filt, RHS8)
    RHS9 = np.dot(Filt, RHS9)


    # Switch to column-vectors (column-major way).
    RHS1col = RHS1.flatten('F')
    RHS2col = RHS2.flatten('F')
    RHS3col = RHS3.flatten('F')
    RHS4col = RHS4.flatten('F')
    RHS5col = RHS5.flatten('F')
    RHS6col = RHS6.flatten('F')
    RHS7col = RHS7.flatten('F')
    RHS8col = RHS8.flatten('F')
    RHS9col = RHS9.flatten('F')

    # hack tha Dirichlet boundaries!!
    RHS2col[vmapW] = 0
    RHS3col[vmapW] = 0

    RHS5col[vmapW] = 0
    RHS6col[vmapW] = 0

    RHS8col[vmapW] = 0
    RHS9col[vmapW] = 0

    RHS1 = np.reshape(RHS1col, (ctx.numLocalPoints, K), order='F')
    RHS2 = np.reshape(RHS2col, (ctx.numLocalPoints, K), order='F')
    RHS3 = np.reshape(RHS3col, (ctx.numLocalPoints, K), order='F')
    RHS4 = np.reshape(RHS4col, (ctx.numLocalPoints, K), order='F')
    RHS5 = np.reshape(RHS5col, (ctx.numLocalPoints, K), order='F')
    RHS6 = np.reshape(RHS6col, (ctx.numLocalPoints, K), order='F')
    RHS7 = np.reshape(RHS7col, (ctx.numLocalPoints, K), order='F')
    RHS8 = np.reshape(RHS8col, (ctx.numLocalPoints, K), order='F')
    RHS9 = np.reshape(RHS9col, (ctx.numLocalPoints, K), order='F')

    return (RHS1, RHS2, RHS3, RHS4,
        RHS5, RHS6, RHS7, RHS8, RHS9)

# Main solver:
if __name__ == '__main__':
    g = 9.81
    r = 0.997
    finalTime = 40.0
    t = 0.0

    meshManager = dg.MeshManager()
    meshManager.readMesh('./input/box.msh')

    # Numerical parameters:
    N = 8
    CFL = 0.45

    filtOrder = 8
    filtCutoff = 0.0*N

    nodes = dg.TriangleNodesProvisioner(N, meshManager)
    nodes.buildFilter(filtCutoff, filtOrder)

    outputter = dg.VtkOutputter(nodes)

    ctx = nodes.dgContext()

    x = ctx.x
    y = ctx.y

    Np = ctx.numLocalPoints
    K = ctx.numElements
    eta2 = 0.5*np.exp(-10*(x*x) -10*(y*y), dtype=np.float64 , order='C')
    eta1 = 0*eta2
    eta3 = eta2

    #eta = -1*(x/1500.0)
    u1   = np.zeros([Np, K], dtype=np.float64, order='C')
    v1   = np.zeros([Np, K], dtype=np.float64, order='C')
    H1   = 7.0*np.ones([Np, K], dtype=np.float64, order='C')
    H2 = 2*np.ones([Np, K], dtype=np.float64, order='C')
    H3 = 10.0*np.ones([Np, K], dtype=np.float64, order='C')

    h1 = H1 - eta2 + eta1
    h1u1 = h1*u1
    h1v1 = h1*v1

    h2 = H2 + eta2 - eta3
    eta1 = 0*eta2
    h2u2 = 0*eta1
    h2v2 = 0*eta1

    B = 0.0*eta1
    h3 = H3 + eta3 - B
    h3u3 = 0.0*eta1
    h3v3 = 0.0*eta1

    Bx = 0.0*eta1
    By = 0.0*eta1

    # setup fields dictionary for outputting.
    fields = dict()
    fields["etatop"] = eta1
    fields["etamid"] = -H1 + eta2
    fields["etabot"] = -H1 - H2 + eta3
    fields["bot"] = -H1 -H2 -H3 + B
    outputter.writeFieldsToFiles(fields, 0)

    c = np.sqrt(g*(h1 + h2 + h3))
    #dt =0.125*0.000724295
    dt = 0.5*0.000724295
    #dt = CFL*dx/np.max(abs(c))

    step = 0
    while t < finalTime:
        
        (RHS1,RHS2,RHS3,RHS4,RHS5,
            RHS6,RHS7,RHS8,RHS9) = sw3lrComputeRHS(h1, h1u1, h1v1, h2, h2u2, h2v2, h3, h3u3, h3v3, g, Bx, By, ctx, r)
    
        # predictor
        h1hat  = h1 + 0.5*dt*RHS1
        h1u1hat = h1u1 + 0.5*dt*RHS2
        h1v1hat = h1v1 + 0.5*dt*RHS3

        h2hat  = h2 + 0.5*dt*RHS4
        h2u2hat = h2u2 + 0.5*dt*RHS5
        h2v2hat = h2v2 + 0.5*dt*RHS6

        h3hat  = h3 + 0.5*dt*RHS7
        h3u3hat = h3u3 + 0.5*dt*RHS8
        h3v3hat = h3v3 + 0.5*dt*RHS9

        (RHS1,RHS2,RHS3,RHS4,RHS5,
            RHS6,RHS7,RHS8,RHS9) = sw3lrComputeRHS(h1hat, h1u1hat, h1v1hat, h2hat, h2u2hat, h2v2hat, h3hat, h3u3hat, h3v3hat, g, Bx, By, ctx, r)
    
        # corrector - Update solution
        h1 += dt*RHS1
        h1u1 += dt*RHS2
        h1v1 += dt*RHS3

        h2 += dt*RHS4
        h2u2 += dt*RHS5
        h2v2 += dt*RHS6

        h3 += dt*RHS7
        h3u3 += dt*RHS8
        h3v3 += dt*RHS9

        # layer numerical adjustments
        eta1 = (h1+h2+h3) - (H1+H2+H3-B)
        eta2 = H1 + eta1 - h1
        eta3 = H2 + eta2 - h2

        u1 = h1u1/h1
        u2 = h2u2/h2
        u3 = h3u3/h3

        v1 = h1v1/h1
        v2 = h2v2/h2
        v3 = h3v3/h3

        Filt = ctx.filter
        eta1 = np.dot(ctx.filter, eta1)
        eta2 = np.dot(ctx.filter, eta2)
        eta3 = np.dot(ctx.filter, eta3)

        h1 = H1 + eta1 - eta2
        h2 = H2 + eta2 - eta3
        h3 = H3 + eta3 - B 

        h1u1 = h1*u1
        h2u2 = h2*u2
        h3u3 = h2*u3

        h1v1 = h1*v1
        h2v2 = h2*v2
        h3v3 = h3*v3

        h1 = h1.flatten("F")
        h2 = h2.flatten("F")
        h3 = h3.flatten("F")

        h1u1 = h1u1.flatten("F")
        h1v1 = h1v1.flatten("F")

        h2u2 = h2u2.flatten("F")
        h2v2 = h2v2.flatten("F")

        h3u3 = h3u3.flatten("F")
        h3v3 = h3v3.flatten("F")


        vmapM = ctx.vmapM
        vmapP = ctx.vmapP

        h1[vmapM] = 0.5*(h1[vmapM] + h1[vmapP])
        h1[vmapP] = h1[vmapM]

        h2[vmapM] = 0.5*(h2[vmapM] + h2[vmapP])
        h2[vmapP] = h2[vmapM]

        h3[vmapM] = 0.5*(h3[vmapM] + h3[vmapP])
        h3[vmapP] = h3[vmapM]
        h1u1[vmapM] = 0.5*(h1u1[vmapM] + h1u1[vmapP])
        h1u1[vmapP] = h1u1[vmapM]

        h1v1[vmapM] = 0.5*(h1v1[vmapM] + h1v1[vmapP])
        h1v1[vmapP] = h1v1[vmapM]

        h2u2[vmapM] = 0.5*(h2u2[vmapM] + h2u2[vmapP])
        h2u2[vmapP] = h2u2[vmapM]

        h2v2[vmapM] = 0.5*(h2v2[vmapM] + h2v2[vmapP])
        h2v2[vmapP] = h2v2[vmapM]

        h3u3[vmapM] = 0.5*(h3u3[vmapM] + h3u3[vmapP])
        h3u3[vmapP] = h3u3[vmapM]

        h3v3[vmapM] = 0.5*(h3v3[vmapM] + h3v3[vmapP])
        h3v3[vmapP] = h3v3[vmapM]

        h1 = np.reshape(h1, (ctx.numLocalPoints, K), order='F')
        h2 = np.reshape(h2, (ctx.numLocalPoints, K), order='F')
        h3 = np.reshape(h3, (ctx.numLocalPoints, K), order='F')

        h1u1 = np.reshape(h1u1, (ctx.numLocalPoints, K), order='F')
        h1v1 = np.reshape(h1v1, (ctx.numLocalPoints, K), order='F')

        h2u2 = np.reshape(h2u2, (ctx.numLocalPoints, K), order='F')
        h2v2 = np.reshape(h2v2, (ctx.numLocalPoints, K), order='F')

        h3u3 = np.reshape(h3u3, (ctx.numLocalPoints, K), order='F')
        h3v3 = np.reshape(h3v3, (ctx.numLocalPoints, K), order='F')

        h_max = np.max(np.abs(h1))
        if h_max > 1e8  or np.isnan(h_max):
            raise Exception("A numerical instability has occurred.")

        t += dt
        step += 1

        # outputting
        if step % 5 == 0:
            print('t=' + str(t))
            eta1 = h1+h2+h3-(H1+H2+H3-B)
            eta2 = H1 + eta1 - h1
            eta3 = H2 + eta2 - h2 
            fields["etatop"] = eta1
            fields["etamid"] = -H1 + eta2
            fields["etabot"] = -H1 - H2 + eta3
            fields["bot"] = -H1 -H2 -H3 + B
            outputter.writeFieldsToFiles(fields, step)
