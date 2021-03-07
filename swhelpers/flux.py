def sw2dComputeFluxes(h, hu, hv, hN, g, H):
    #h equation
    F1 = hu
    G1 = hv

    # Get velocity fields
    u = hu / h
    v = hv / h

    # hu equation
    F2 = hu*u + 0.5*g*h*h
    G2 = hu*v

    # hv equation
    F3 = hv*u
    G3 = hv*v + 0.5*g*h*h

    # N (tracer) equation
    F4 = hN*u
    G4 = hN*v

    return ((F1,F2,F3,F4),(G1,G2,G3,G4))