import numpy as np

def positivityPreservingLimiter2D(h, hu, hv):
    Np, K = h.shape
    hmin = np.tile(np.min(h, axis=0), (Np, 1))
    hmin[hmin < 1e-3] = 1e-3

    hmean = np.tile(np.mean(h, axis=0), (Np, 1))

    theta = np.ones((Np, K))
    theta = hmean / (hmean - hmin + 1e-4)
    
    theta[theta > 1] = 1.0
    humean = np.tile(np.mean(hu, axis=0), (Np, 1))
    hvmean = np.tile(np.mean(hv, axis=0), (Np, 1))

    h  = hmean  + theta*(h  - hmean)
    hu = humean + theta*(hu - humean)
    hv = hvmean + theta*(hv - hvmean)

    return h, hu, hv


def minmod(a, b):
    soln = np.zeros(a.shape)
    for i, _ in enumerate(a):
        if a[i] < b[i] and a[i]*b[i] > 0:
            soln[i] = a[i]
        elif b[i] < a[i] and a[i]*b[i] > 0:
            soln[i] = b[i]
        else:
            soln[i] = 0.0
    
    return soln

def surfaceReconstruction(etaM, hM, etaP, hP):
    # get bed elevations
    zM = etaM - hM
    zP = etaP - hP

    dz = (zP - 0.5*minmod(zP - zM, 1e-3*np.ones(zM.shape))) - (zM + 0.5*minmod(zM-zP, -1e-3*np.ones(zM.shape)))


    # flux limit
    #etaCorrM = zP - zM - dz
    #for i,_ in enumerate(etaCorrM):
    #    if etaCorrM[i] > (etaP[i] - etaM[i]):
    #        etaCorrM[i] = etaP[i] - etaM[i]

    #    if etaCorrM[i] < 0:
    #        etaCorrM[0] = 0.0    
    
    #etaM += etaCorrM

    etaCorrP = zM - zP - dz
    for i, _ in enumerate(etaCorrP):
        if etaCorrP[i] > (etaM[i] - etaP[i]):
            etaCorrP[i] = etaM[i] - etaP[i]
            
        if etaCorrP[i] <= 0:
            etaCorrP[i] = 0.0
        else:
            etaP[i] += etaCorrP[i]


    # Get corrected bed elevation
    #zM = etaM - hM
    zP = etaP - hP

    # enforce non-negativity
    maxz = zM
    for i, _ in enumerate(zM):
        if zP[i] > zM[i]:
            maxz[i] = zP[i]

    hM = etaM - maxz
    hM[hM <= 1e-3] = 1e-3
    hP = etaP - maxz
    hP[hP <= 1e-3] = 1e-3

    return hM, hP