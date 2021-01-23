import numpy as np
from scipy.linalg import solve_triangular

from .flux import sw2dComputeFluxes

def sw2dComputeRHS_curved(h, hu, hv, hN, zx, zy, g, H, f, CD, ctx, cub_ctx, gauss_ctx, curvedEls, J, gmapM, gmapP):

    cub_h  = np.dot(cub_ctx.V, h)
    cub_hu = np.dot(cub_ctx.V, hu)
    cub_hv = np.dot(cub_ctx.V, hv)
    cub_hN = np.dot(cub_ctx.V, hN)

    cub_zx = np.dot(cub_ctx.V, zx)
    cub_zy = np.dot(cub_ctx.V, zy)

    cub_H = np.dot(cub_ctx.V, H)

    ((F1,F2,F3,F4),(G1,G2,G3,G4)) = sw2dComputeFluxes(cub_h, cub_hu, cub_hv, cub_hN, g, cub_H)

    DrT = np.transpose(cub_ctx.Dr) 
    DsT = np.transpose(cub_ctx.Ds)

    tmpr = cub_ctx.W*(cub_ctx.rx*F1 + cub_ctx.ry*G1)
    tmps = cub_ctx.W*(cub_ctx.sx*F1 + cub_ctx.sy*G1)
    ddr = np.dot(DrT, tmpr)
    dds = np.dot(DsT, tmps)
    MMRHS1 = ddr + dds

    tmpr = cub_ctx.W*(cub_ctx.rx*F2 + cub_ctx.ry*G2)
    tmps = cub_ctx.W*(cub_ctx.sx*F2 + cub_ctx.sy*G2)
    ddr = np.dot(DrT, tmpr)
    dds = np.dot(DsT, tmps)
    MMRHS2 = ddr + dds

    tmpr = cub_ctx.W*(cub_ctx.rx*F3 + cub_ctx.ry*G3)
    tmps = cub_ctx.W*(cub_ctx.sx*F3 + cub_ctx.sy*G3)
    ddr = np.dot(DrT, tmpr)
    dds = np.dot(DsT, tmps)
    MMRHS3 = ddr + dds

    tmpr = cub_ctx.W*(cub_ctx.rx*F4 + cub_ctx.ry*G4)
    tmps = cub_ctx.W*(cub_ctx.sx*F4 + cub_ctx.sy*G4)
    ddr = np.dot(DrT, tmpr)
    dds = np.dot(DsT, tmps)
    MMRHS4 = ddr + dds

    nx = gauss_ctx.nx
    ny = gauss_ctx.ny
    mapW = gauss_ctx.BCmap[3]

    nxW = nx.flatten('F')[mapW]
    nyW = ny.flatten('F')[mapW]

    NGauss      = gauss_ctx.nx.shape[0]
    numElements = gauss_ctx.nx.shape[1]

    gauss_h = np.dot(gauss_ctx.Interp, h).flatten('F')
    hM = gauss_h[gmapM]
    hP = gauss_h[gmapP]

    gauss_hu = np.dot(gauss_ctx.Interp, hu).flatten('F')
    huM = gauss_hu[gmapM]
    huP = gauss_hu[gmapP]

    gauss_hv = np.dot(gauss_ctx.Interp, hv).flatten('F')
    hvM = gauss_hv[gmapM]
    hvP = gauss_hv[gmapP]

    gauss_hN = np.dot(gauss_ctx.Interp, hN).flatten('F')
    hNM = gauss_hN[gmapM]
    hNP = gauss_hN[gmapP]

    gauss_H = np.dot(gauss_ctx.Interp, H).flatten('F')
    HM = gauss_H[gmapM]
    HP = gauss_H[gmapP]

    uM = huM / hM
    uP = huP / hP

    vM = hvM / hM
    vP = hvP / hP

    huP[mapW] = huM[mapW] - 2*nxW*(huM[mapW]*nxW + hvM[mapW]*nyW)
    hvP[mapW] = hvM[mapW] - 2*nyW*(huM[mapW]*nxW + hvM[mapW]*nyW)

    ((F1M,F2M,F3M,F4M),(G1M,G2M,G3M,G4M)) = sw2dComputeFluxes(hM, huM, hvM, hNM, g, HM)
    ((F1P,F2P,F3P,F4P),(G1P,G2P,G3P,G4P)) = sw2dComputeFluxes(hP, huP, hvP, hNP, g, HP)

    spdM = np.sqrt(uM*uM + vM*vM) + np.sqrt(g*hM)
    spdP = np.sqrt(uP*uP + vP*vP) + np.sqrt(g*hP)

    spdMax = np.max(np.array([spdM, spdP]), axis=0)

    numFaces = 3
    Nfp = int(NGauss / numFaces)
    K = ctx.numElements
    
    lam = np.reshape(spdMax, (Nfp, numFaces*numElements), order='F')
    lamMaxMat = np.outer(np.ones((Nfp, 1), dtype=np.float), np.max(lam, axis=0))
    spdMax = np.reshape(lamMaxMat, (Nfp*numFaces, K), order='F')

    dh = np.reshape(hM - hP, (Nfp*numFaces, K), order='F')
    dhu = np.reshape(huM - huP, (Nfp*numFaces, K), order='F')
    dhv = np.reshape(hvM - hvP, (Nfp*numFaces, K), order='F')
    dhN = np.reshape(hNM - hNP, (Nfp*numFaces, K), order='F')

    F1M = np.reshape(F1M, (Nfp*numFaces, K), order='F')
    F2M = np.reshape(F2M, (Nfp*numFaces, K), order='F')
    F3M = np.reshape(F3M, (Nfp*numFaces, K), order='F')
    F4M = np.reshape(F4M, (Nfp*numFaces, K), order='F')

    G1M = np.reshape(G1M, (Nfp*numFaces, K), order='F')
    G2M = np.reshape(G2M, (Nfp*numFaces, K), order='F')
    G3M = np.reshape(G3M, (Nfp*numFaces, K), order='F')
    G4M = np.reshape(G4M, (Nfp*numFaces, K), order='F')

    F1P = np.reshape(F1P, (Nfp*numFaces, K), order='F')
    F2P = np.reshape(F2P, (Nfp*numFaces, K), order='F')
    F3P = np.reshape(F3P, (Nfp*numFaces, K), order='F')
    F4P = np.reshape(F4P, (Nfp*numFaces, K), order='F')

    G1P = np.reshape(G1P, (Nfp*numFaces, K), order='F')
    G2P = np.reshape(G2P, (Nfp*numFaces, K), order='F')
    G3P = np.reshape(G3P, (Nfp*numFaces, K), order='F')
    G4P = np.reshape(G4P, (Nfp*numFaces, K), order='F')

    fluxh  = 0.5*((F1M + F1P)*nx + (G1M+G1P)*ny + spdMax*dh)
    fluxhu = 0.5*((F2M + F2P)*nx + (G2M+G2P)*ny + spdMax*dhu)
    fluxhv = 0.5*((F3M + F3P)*nx + (G3M+G3P)*ny + spdMax*dhv)
    fluxhN = 0.5*((F4M + F4P)*nx + (G4M+G4P)*ny + spdMax*dhN)

    interp = np.transpose(gauss_ctx.Interp)
    MMRHS1 -= np.dot(interp, (gauss_ctx.W*fluxh))   # this was weird...
    MMRHS2 -= np.dot(interp, (gauss_ctx.W*fluxhu))
    MMRHS3 -= np.dot(interp, (gauss_ctx.W*fluxhv))
    MMRHS4 -= np.dot(interp, (gauss_ctx.W*fluxhN))

    RHS1 = np.zeros((ctx.numLocalPoints, ctx.numElements))
    RHS2 = np.zeros((ctx.numLocalPoints, ctx.numElements))
    RHS3 = np.zeros((ctx.numLocalPoints, ctx.numElements))
    RHS4 = np.zeros((ctx.numLocalPoints, ctx.numElements))

    curvedElsSet = set(curvedEls)
    straightEls = []
    for k in range(0, numElements):
        if k not in curvedElsSet:
            straightEls.append(k)
    
    V = ctx.V
    VT = V.T
    mmInvStandard = np.dot(V, VT)
    RHS1 = np.dot(mmInvStandard, MMRHS1/J)
    RHS2 = np.dot(mmInvStandard, MMRHS2/J)
    RHS3 = np.dot(mmInvStandard, MMRHS3/J)
    RHS4 = np.dot(mmInvStandard, MMRHS4/J)

    for k in curvedElsSet:
        mmChol = cub_ctx.MMChol[:,:,k]
        RHS1[:, k] = solve_triangular(mmChol, solve_triangular(mmChol, MMRHS1[:, k], trans='T'))
        RHS2[:, k] = solve_triangular(mmChol, solve_triangular(mmChol, MMRHS2[:, k], trans='T'))
        RHS3[:, k] = solve_triangular(mmChol, solve_triangular(mmChol, MMRHS3[:, k], trans='T'))
        RHS4[:, k] = solve_triangular(mmChol, solve_triangular(mmChol, MMRHS4[:, k], trans='T'))

    # Add source terms
    u = hu / h
    v = hv / h

    norm_u = np.hypot(u, v)
    CD_norm_u = CD*norm_u

    RHS2 += f*hv - CD_norm_u*u
    RHS3 -= f*hu - CD_norm_u*v
    RHS2 -= g*h*zx   # note: should over-integrate these.
    RHS3 -= g*h*zy

    return (RHS1, RHS2, RHS3, RHS4)

def sw2dComputeRHS(h, hu, hv, hN, zx, zy, g, H, f, CD, ctx, vmapM, vmapP):
    BCmap = ctx.BCmap
    nx = ctx.nx
    ny = ctx.ny
    rx = ctx.rx
    sx = ctx.sx
    ry = ctx.ry
    sy = ctx.sy
    Dr = ctx.Dr
    Ds = ctx.Ds
    Nfp = ctx.numFacePoints
    K = ctx.numElements

    Lift = ctx.Lift
    Fscale = ctx.Fscale

    hC = h.flatten('F')
    huC = hu.flatten('F')
    hvC = hv.flatten('F')
    hNC = hN.flatten('F')
    nxC = nx.flatten('F')
    nyC = ny.flatten('F')
    
    mapW = BCmap[3]

    # get field values along elemental faces.
    hM = hC[vmapM]
    hP = hC[vmapP]

    eta = h - H
    etaC = eta.flatten('F')
    #etaM = etaC[vmapM]
    #etaP = etaC[vmapP]

    uM = huC[vmapM] / hC[vmapM]
    uP = huC[vmapP] / hC[vmapP]

    vM = hvC[vmapM] / hC[vmapM]
    vP = hvC[vmapP] / hC[vmapP]

    hNM = hNC[vmapM]
    hNP = hNC[vmapP]

    nxW = nxC[mapW]
    nyW = nyC[mapW]

    # hM, hP = surfaceReconstruction(etaM, hM, etaP, hP)
    # h = np.reshape(hC, (Np, K), order='F')

    # re-form conserved transport from corrected 
    # water column heights.
    huM = hM*uM
    hvM = hM*vM

    huP = hP*uP
    hvP = hP*vP

    # set bc's (no normal flow thru the walls).
    huP[mapW] = huM[mapW] - 2*nxW*(huM[mapW]*nxW + hvM[mapW]*nyW)
    hvP[mapW] = hvM[mapW] - 2*nyW*(huM[mapW]*nxW + hvM[mapW]*nyW)

    # compute jump in states
    dh = hM - hP
    dhu = huM - huP
    dhv = hvM - hvP
    dhN = hNM - hNP

    ((F1M,F2M,F3M,F4M),(G1M,G2M,G3M,G4M)) = sw2dComputeFluxes(hM, huM, hvM, hNM, g, H)
    ((F1P,F2P,F3P,F4P),(G1P,G2P,G3P,G4P)) = sw2dComputeFluxes(hP, huP, hvP, hNP, g, H)
    ((F1,F2,F3,F4),(G1,G2,G3,G4)) = sw2dComputeFluxes(h, hu, hv, hN, g, H)

    uM = huM/hM 
    vM = hvM/hM

    uP = huP/hP
    vP = hvP/hP

    spdM = np.sqrt(uM*uM + vM*vM) + np.sqrt(g*hM)
    spdP = np.sqrt(uP*uP + vP*vP) + np.sqrt(g*hP)

    spdMax = np.max(np.array([spdM, spdP]), axis=0)

    # spdMax = np.max(spdMax)
    lam = np.reshape(spdMax, (ctx.numFacePoints, ctx.numFaces*ctx.numElements), order='F')
    lamMaxMat = np.outer(np.ones((Nfp, 1), dtype=np.float), np.max(lam, axis=0))
    spdMax = lamMaxMat.flatten('F')

    # strong form: Compute flux jump vector. (fluxM - numericalFlux ) dot nW
    dFlux1 = 0.5*((F1M - F1P)*nxC + (G1M-G1P)*nyC - spdMax*dh)
    dFlux2 = 0.5*((F2M - F2P)*nxC + (G2M-G2P)*nyC - spdMax*dhu)
    dFlux3 = 0.5*((F3M - F3P)*nxC + (G3M-G3P)*nyC - spdMax*dhv)
    dFlux4 = 0.5*((F4M - F4P)*nxC + (G4M-G4P)*nyC - spdMax*dhN)

    dFlux1Mat = np.reshape(dFlux1, (Nfp*ctx.numFaces, K), order='F')
    dFlux2Mat = np.reshape(dFlux2, (Nfp*ctx.numFaces, K), order='F')
    dFlux3Mat = np.reshape(dFlux3, (Nfp*ctx.numFaces, K), order='F')
    dFlux4Mat = np.reshape(dFlux4, (Nfp*ctx.numFaces, K), order='F')

    # Flux divergence:
    RHS1 = -(rx*np.dot(Dr, F1) + sx*np.dot(Ds, F1))
    RHS1+= -(ry*np.dot(Dr, G1) + sy*np.dot(Ds, G1))

    RHS2 = -(rx*np.dot(Dr, F2) + sx*np.dot(Ds, F2))
    RHS2+= -(ry*np.dot(Dr, G2) + sy*np.dot(Ds, G2))

    RHS3 = -(rx*np.dot(Dr, F3) + sx*np.dot(Ds, F3))
    RHS3+= -(ry*np.dot(Dr, G3) + sy*np.dot(Ds, G3))

    RHS4 = -(rx*np.dot(Dr, F4) + sx*np.dot(Ds, F4))
    RHS4+= -(ry*np.dot(Dr, G4) + sy*np.dot(Ds, G4))

    surfaceRHS1 = Fscale*dFlux1Mat
    surfaceRHS2 = Fscale*dFlux2Mat
    surfaceRHS3 = Fscale*dFlux3Mat
    surfaceRHS4 = Fscale*dFlux4Mat

    RHS1 += np.dot(Lift, surfaceRHS1)
    RHS2 += np.dot(Lift, surfaceRHS2)
    RHS3 += np.dot(Lift, surfaceRHS3)
    RHS4 += np.dot(Lift, surfaceRHS4)

    # Add source terms
    u = hu / h
    v = hv / h

    norm_u = np.hypot(u, v)
    CD_norm_u = CD*norm_u

    RHS2 += f*hv - CD_norm_u*u
    RHS3 -= f*hu - CD_norm_u*v
    RHS2 -= g*h*zx
    RHS3 -= g*h*zy

    return (RHS1, RHS2, RHS3, RHS4)