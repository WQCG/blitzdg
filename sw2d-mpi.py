#!/usr/bin/python3
# Run with: mpirun -np <#> python3 sw2d-mpi.py
'''
Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
See COPYING and LICENSE files at project root for more details.
'''

import numpy as np
from pyblitzdg import pyblitzdg as dg
from mpi4py import MPI as mpi
import sw2d as sw2d
from matplotlib import pyplot as plt

def master():
    return mpi.COMM_WORLD.Get_rank() == 0

def sw2dmpiComputeRHS(h, hu, hv, g, H, ctx, commBdryNodes):
    vmapM = ctx.vmapM
    vmapP = ctx.vmapP
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
    
    mapW = ctx.BCmap[3]

    # get field values along elemental faces.
    hM = hC[vmapM]
    hP = hC[vmapP]

    huM = huC[vmapM]
    huP = huC[vmapP]

    hvM = hvC[vmapM]
    hvP = hvC[vmapP]

    nxW = nxC[mapW]
    nyW = nyC[mapW]

    # Rest of code is same as if serial...
    # set bc's (no normal flow thru the walls).
    huP[mapW] = huM[mapW] - 2*nxW*(huM[mapW]*nxW + hvM[mapW]*nyW)
    hvP[mapW] = hvM[mapW] - 2*nyW*(huM[mapW]*nxW + hvM[mapW]*nyW)


    xM = ctx.x.flatten('F')[vmapM]
    yM = ctx.y.flatten('F')[vmapM]

    xW = xM[mapW]
    yW = yM[mapW]
    hW = hM[mapW]
    huW = huM[mapW]
    hvW = hvM[mapW]



    # populate ghost nodes with neighboring sub-domain's values
    comm = mpi.COMM_WORLD
    numProcs = comm.Get_size()
    rank = comm.Get_rank()
    sends  = []
    recvs = []
    recvsData = [None]*numProcs
    for i in [r for r in range(0, numProcs) if r != rank]:
             
        if len(commBdryNodes[i]) > 0:
            commData = np.ndarray((len(commBdryNodes[i]), 5), dtype=np.dtype('Float64') )
            commData[:, 0] = hW[commBdryNodes[i]]
            commData[:, 1] = huW[commBdryNodes[i]]
            commData[:, 2] = hvW[commBdryNodes[i]]
            commData[:, 3] = xW[commBdryNodes[i]]
            commData[:, 4] = yW[commBdryNodes[i]]

            print(f"Rank {rank}: async send to {i}")
            req = comm.Isend([commData, mpi.DOUBLE], i, rank*100 +i )
            sends.append(req)

        recvsData[i] = np.ndarray((len(commBdryNodes[i]), 5), dtype=np.dtype('Float64'))
        recvs.append(comm.Irecv([recvsData[i], mpi.DOUBLE], source=i,  tag=i*100+rank))

    if len(sends) > 0:
        mpi.Request.waitall(sends)

    # could speed this up with a 'wait any' and unpack, until received all.
    if len(recvs) > 0:
        mpi.Request.waitall(recvs)

    #print(f"[{rank}]: {ctx.x.flatten('F')[vmapM][mapComm]}")
    #exit(0)
    
    
    # Populate 'internal' ghost states with received data  
    for i in [r for r in range(0, numProcs) if r != rank]:

        mapComm = commBdryNodes[i]
        if len(mapComm) > 0:
            xCommRecv = recvsData[i][:, 3]
            yCommRecv = recvsData[i][:, 4]
            
            pointsRecv = np.c_[xCommRecv, yCommRecv]

            xCommSend = xW[mapComm]
            yCommSend = yW[mapComm]

            pointsSend = np.c_[xCommSend, yCommSend]

            #for ind in range(0,len(xCommRecv)):
            #    print(f"[{rank}]: ({xCommRecv[ind]}, {yCommRecv[ind]})")

            #if not master():
                #plt.figure()
                #plt.hold(True)
                #for z in range(0, len(xCommRecv)-1):
                    #plt.plot(xCommRecv[z:z+2], yCommRecv[z:z+2], 'r') # , xCommSend, yCommSend, '.g'
                    #plt.axis([-1, 1, -1, 1])
                    #plt.draw()
                    #plt.show()
                    #plt.ioff()
                    #plt.pause(0.25)
                    #plt.ioff()

            
            # this seems to work, then...


            commRecvMap = np.ndarray((len(xCommSend), 2), dtype=np.dtype('Int32'))
            for j in range(0, len(xCommSend)):
                dists1 = np.hypot(xCommSend[j] - xCommRecv, yCommSend[j] - yCommRecv)
                #dists2 = np.hypot(xCommSend[j+1] - xCommRecv, yCommSend[j+1] - yCommRecv) 
                mask1 = dists1 < 1e-5
                #mask2 = dists2 < 1e-5
                inds1 = np.arange(0, len(xCommRecv))[mask1]
                #inds2 = np.arange(0, len(xCommRecv))[mask2]

                #if np.hypot(xCommSend[j] - xCommSend[j+1], yCommSend[j] - yCommSend[j+1]) -
                #    np.hypot(xCommRecv[inds1[0]]) 
                
                # DEREK: fix this later
                commRecvMap[j, 0] = inds1[0]
                commRecvMap[j, 1] = inds1[1]

            print(f"[{rank}]: {commRecvMap}")
        #print(f"[{rank}]: {pointsSend}")
        #exit(0)
        

   
    for i in [r for r in range(0, numProcs) if r != rank]:
        mapComm = commBdryNodes[i]
        # Derekian diffusion, should be fixed.
        #hP[mapW][mapComm]  = np.min(recvsData[i][commRecvMap[:, 0].flatten(), 0], recvsData[i][commRecvMap[:, 1].flatten(), 0])
        #huP[mapW][mapComm] = np.min(recvsData[i][commRecvMap[:, 0].flatten(), 1], recvsData[i][commRecvMap[:, 1].flatten(), 1])
        #hvP[mapW][mapComm] = np.min(recvsData[i][commRecvMap[:, 0].flatten(), 2], recvsData[i][commRecvMap[:, 1].flatten(), 2])
            
    # compute jump in states
    dh = hM - hP
    dhu = huM - huP
    dhv = hvM - hvP

    ((F1M,F2M,F3M),(G1M,G2M,G3M)) = sw2d.sw2dComputeFluxes(hM, huM, hvM, g, H)
    ((F1P,F2P,F3P),(G1P,G2P,G3P)) = sw2d.sw2dComputeFluxes(hP, huP, hvP, g, H)
    ((F1,F2,F3),(G1,G2,G3)) = sw2d.sw2dComputeFluxes(h, hu, hv, g, H)

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

    K = ctx.numElements
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

comm = mpi.COMM_WORLD
numProcs = comm.Get_size()
rank = comm.Get_rank()
name = mpi.Get_processor_name()

# Main solver:
g = 9.81
finalTime = 40.0
t = 0.0

# all procs need to do this...
meshManager = dg.MeshManager()
meshManager.readMesh('./input/box.msh')
vertices = meshManager.vertices
elements = meshManager.elements

Nv = len(vertices)
Ne = len(elements)

# default vertex/elements maps -- all zeros (no communications).
vertexMap = np.zeros(Nv, dtype=np.dtype("Int32"))
elementMap = np.zeros(Ne, dtype=np.dtype("Int32"))

# partition the mesh on the master node
if master():
    if numProcs >= 2:
        meshManager.partitionMesh(numProcs)

        vertexMap = meshManager.vertexPartitionMap
        elementMap = meshManager.elementPartitionMap

        for p in range(1, numProcs):
            comm.Send([vertexMap, mpi.INT], dest=p, tag=10+p)
            comm.Send([elementMap, mpi.INT], dest=p, tag=4096+p)

elif not master() and numProcs >= 2:
    vertexMap = np.empty(Nv, dtype=np.dtype('Int32'))
    elementMap = np.empty(Ne, dtype=np.dtype('Int32'))
    comm.Recv([vertexMap, mpi.INT], source=0, tag=10+rank)
    comm.Recv([elementMap, mpi.INT], source=0, tag=4096+rank)


print(vertexMap.shape)
print(vertices.shape)
vertsLocal        = vertices[vertexMap==rank,  :]
elemsLocal2Global = elements[elementMap==rank, :]

# build map from local to global vertex numbers
vertIndLocal2Global = np.arange(0, Nv)[vertexMap==rank]
elemIndLocal2Global = np.arange(0, Ne)[elementMap==rank]

# build reverse-lookups tables as well.
vertIndGlobal2Local = dict.fromkeys(vertIndLocal2Global)
for i, k in enumerate(vertIndGlobal2Local.keys()):
    vertIndGlobal2Local[k] = i

elemIndGlobal2Local = dict.fromkeys(elemIndLocal2Global)
for i, k in enumerate(elemIndGlobal2Local.keys()):
    elemIndGlobal2Local[k] = i

# get number of faces per element
numFaces = elements.shape[1]

# Set up local EToV
NeLocal = len(elemIndLocal2Global)
NvLocal = len(vertIndLocal2Global)
EToV = np.zeros((NeLocal, numFaces))

commVerts = np.zeros((NeLocal, numFaces)) - np.Infinity
for e in range(0, NeLocal):
    badVerts = []
    for v in range(0, numFaces):
        eG = elemIndLocal2Global[e]
        vG = elements[eG, v]
        if vG not in vertIndGlobal2Local.keys():
            print(f"[WARNING]: {vG} not found in reverse look-ups. Rank: {rank}. Adding it to this subdomain." )

            # append new row.
            vertsLocal = np.r_[vertsLocal, [vertices[vG, :]]]
            vertIndGlobal2Local[vG] = len(vertsLocal) - 1
            
            badVerts.append((eG, v))
            print(f"Rank: {rank}. badVerts: {badVerts}")

            # Get the subdomain number of the vertex we didn't know about,
            # that will tell us which subdomain to talk to.
            commVerts[e, v] = vertexMap[vG]
            
        vL = vertIndGlobal2Local[vG]
        EToV[e, v] = vL

# de-allocate the global guy.
meshManager = None

localMeshManager = dg.MeshManager()

print("building mesh")
localMeshManager.buildMesh(EToV, vertsLocal)


print("building nodes")
nodes = dg.TriangleNodesProvisioner(1, localMeshManager)
ctx = nodes.dgContext()

BCmap = ctx.BCmap
mapW = BCmap[3]

# check if our wall coincides with our neighbor's wall, if it does, then it's an interior boundary.

xW = ctx.x.flatten('F')[ctx.vmapM][mapW]
yW = ctx.y.flatten('F')[ctx.vmapM][mapW]

sends = []
recvs = []
recvsData = [None]*numProcs

for p in [r for r in range(0, numProcs) if r != rank]:
    commData = np.ndarray((len(mapW), 2), dtype=np.dtype('Float64') )
    commData[:, 0] = xW
    commData[:, 1] = yW
    
    print(f"Rank {rank}: async send to {p}")
    req = comm.Isend([commData, mpi.DOUBLE], p, tag=rank*100 +p )
    sends.append(req)

mpi.Request.waitall(sends)

for p in [r for r in range(0, numProcs) if r != rank]:
    recvsData[p] = np.ndarray((len(mapW), 2), dtype=np.dtype('Float64'))
    recvs.append(comm.Irecv([recvsData[p], mpi.DOUBLE], source=p,  tag=p*100+rank))
    print(f"Rank {rank}: async recv from {p}")


print(recvs)
mpi.Request.waitall(recvs)

points = np.c_[xW, yW]

commBdryNodes = dict()
for p in [r for r in range(0, numProcs) if r != rank]:
    xWrecv = recvsData[p][:, 0]
    yWrecv = recvsData[p][:, 1]

    pointsRecv = np.c_[xWrecv, yWrecv]

    ovLap = []
    for i in np.arange(0, len(xW)):
        dists = np.hypot(points[i, 0] - pointsRecv[:, 0], points[i, 1] - pointsRecv[:, 1]) 
        if np.any([dists <= 1e-9]):
            ovLap.append(i)

    commBdryNodes[p] = ovLap

#if not master():
#    plt.figure()
#    plt.plot(xW[ovLap], yW[ovLap]) # 'og', xWrecv, yWrecv, '.r'
#    plt.draw()
#    plt.show()

print("building outputter")
outputter = dg.VtkOutputter(nodes)
x = ctx.x
y = ctx.y

Np = ctx.numLocalPoints
K = ctx.numElements

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
fields["eta_" + str(rank) + "_"] = eta
fields["u_" + str(rank) + "_"] = u
fields["v_" + str(rank) + "_"] = v
outputter.writeFieldsToFiles(fields, 0)

c = np.sqrt(g*h)
#dt =0.000724295
CFL = 0.1
dx = 0.01
dt = CFL*dx/np.max(abs(c))

step = 0
while t < finalTime:

    (RHS1,RHS2,RHS3) = sw2dmpiComputeRHS(h, hu, hv, g, H, ctx, commBdryNodes)
    
    # predictor
    h1  = h + 0.5*dt*RHS1
    hu1 = hu + 0.5*dt*RHS2
    hv1 = hv + 0.5*dt*RHS3

    (RHS1,RHS2,RHS3) = sw2dmpiComputeRHS(h1, hu1, hv1, g, H, ctx, commBdryNodes)

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
    fields["eta_" + str(rank) + "_"] = eta
    fields["u_" + str(rank) + "_"] = hu/h
    fields["v_" + str(rank) + "_"] = hv/h
    outputter.writeFieldsToFiles(fields, step)
