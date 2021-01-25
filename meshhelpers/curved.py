import numpy as np
import pyblitzdg as dg
from scipy.interpolate import splev, splrep, interp1d, griddata

def adjustStraightEdges(Verts, EToV, bcFaces, xSpline, ySpline, ctx):
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

        hyps1 = np.hypot(xSpline - v1[0], ySpline-v1[1])
        hyps2 = np.hypot(xSpline - v2[0], ySpline-v2[1])

        minInd1 = np.argmin(hyps1)
        minInd2 = np.argmin(hyps2)

        mytol = 200

        if hyps1[minInd1] > mytol or hyps2[minInd2] > mytol:
            # nothing to do here
            continue

        curvedFaces.append([el, face])

        # set new vertex coordinates to closest spline point coordinates
        newx1 = xSpline[minInd1] 
        newy1 = ySpline[minInd1]

        newx2 = xSpline[minInd2]
        newy2 = ySpline[minInd2]

        # update mesh vertex locations
        Verts[v1ind, 0] = newx1
        Verts[v1ind, 1] = newy1

        Verts[v2ind, 0] = newx2
        Verts[v2ind, 1] = newy2

        modifiedVerts[v1ind] = 1  
        modifiedVerts[v2ind] = 1

    return Verts, modifiedVerts, curvedFaces

def deformAndBlendElements(Verts, EToV, curvedFaces, xSpline, ySpline, ss, xs, ys, ctx, NOrder):
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

        v1_dists2 = (x1-xSpline)**2 + (y1-ySpline)**2
        v2_dists2 = (x2-xSpline)**2 + (y2-ySpline)**2

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
        splxev = splev(tLGL, xs, ext=3)
        splyev = splev(tLGL, ys, ext=3) 
        x1 = ctx.x[fmsk, k]
        y1 = ctx.y[fmsk, k]
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

        ctx.x[ids, k] += blend*vdx[ids]
        ctx.y[ids, k] += blend*vdy[ids]
    
    return ctx.x, ctx.y, curvedEls