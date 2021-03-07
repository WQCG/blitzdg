import numpy as np

def makeMapsPeriodic(vmapM, vmapP, vmapO, xFlat, yFlat, xO, yO):
    vmapOM = []
    vmapOP = []
    for i in vmapO:
        xi = xFlat[i]
        yi = yFlat[i]

        ips = np.where( np.logical_and(np.abs(yO - yi) < 1e-3, np.abs(xO - xi) > 1000))[0]

        vmapOM.append(i)
        vmapOP.append(vmapO[ips])

    lookup = dict.fromkeys(vmapOM)

    # Figure out look-up for periodic BCs
    vmapOPflat = []
    for i, l in enumerate(vmapOP):
        if len(l) == 1:
            vmapOPflat.append(l[0])
        else:
            if i==0 and (abs(vmapOP[i+1] -l[0]) <= 1):
                vmapOPflat.append(l[0])
            elif (i > 0 and i < len(vmapOP)-1) and (abs(vmapOP[i-1][0] -l[0]) <= 1 or abs(vmapOP[i+1][0] -l[0]) <= 1):
                vmapOPflat.append(l[0])
            elif (i==len(vmapOP)-1) and (abs(vmapOP[i-1] -l[0]) <= 1):
                vmapOPflat.append(l[0])
            else:
                vmapOPflat.append(l[1])

    vmapOP = vmapOPflat


    for i, key in enumerate(lookup.keys()):
        lookup[key] = vmapOP[i]

    # Need to mutate vmapP
    for i  in range(0, len(vmapM)):
        if vmapM[i] in lookup.keys():
            vmapP[i] = lookup[vmapM[i]]

    return vmapM, vmapP

# Correct boundary condition table for the different types
# TODO, store magic integers like 6 - Dirichlet and 7 - Neuman in an enum prop.
def correctBCTable(bcType, EToV, Verts, bcTag):
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

            if  np.abs(midx - 0.0) < 1e-6:
                bcType[e, faceInd] = bcTag   # outflow
            elif np.abs(midx - 8000.0) < 1e-6:
                bcType[e, faceInd] = bcTag   # outflow

    return bcType