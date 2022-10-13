
function computeThirdBodyPerterbations(cart, t, p::AbstractEOMParams)
    # Compute shifted time
    ts  = p.initEpoch + t

    # Compute third body accelerations
    ax  = 0.0
    ay  = 0.0
    az  = 0.0
    for i in eachindex(p.thirdBodyEphemerides.targIDs)
        # Grab current ephemeride
        ephem   = p.thirdBodyEphemerides.ephems[i]

        # Get gravitational parameter
        μ   = AstroUtils.getTargetGM(ephem)

        # Get position
        s   = getPosition(ephem, ts)

        # Compute relative position of spacecraft with respect to the third body
        d   = SVector(cart[1] - s[1], cart[2] - s[2], cart[3] - s[3])

        # Compute Battin's F(q) function
        Fqk = computeFk(view(cart, 1:3), s)

        # Compute acceleration
        dk3     = sqrt(d[1]*d[1] + d[2]*d[2] + d[3]*d[3])^3
        dk3Inv  = 1.0 / dk3
        ax     -= μ*dk3Inv*(cart[1] + Fqk*s[1])
        ay     -= μ*dk3Inv*(cart[2] + Fqk*s[2])
        az     -= μ*dk3Inv*(cart[3] + Fqk*s[3])
    end

    # Currently does nothing
    return SVector(ax, ay, az)
end

# Computes the third body perterbations and partials w.r.t cartesian state (add time partials in future if needed)
function computeThirdBodyPerterbationsAndPartials(cartWithMass, t, p::AbstractEOMParams)
    # Compute shifted time
    ts  = p.initEpoch + t

    # Get individual componants of full state vector
    cart    = view(cartWithMass, 1:6)
    m       = cartWithMass[7]

    # Compute third body accelerations
    ax      = 0.0
    ay      = 0.0
    az      = 0.0
    daxdrx  = 0.0
    daxdry  = 0.0
    daxdrz  = 0.0
    daydrx  = 0.0
    daydry  = 0.0
    daydrz  = 0.0
    dazdrx  = 0.0
    dazdry  = 0.0
    dazdrz  = 0.0
    daxdvx  = 0.0
    daxdvy  = 0.0
    daxdvz  = 0.0
    daydvx  = 0.0
    daydvy  = 0.0
    daydvz  = 0.0
    dazdvx  = 0.0
    dazdvy  = 0.0
    dazdvz  = 0.0
    daxdm   = 0.0
    daydm   = 0.0
    dazdm   = 0.0
    for i in eachindex(p.thirdBodyEphemerides.targIDs)
        # Grab current ephemeride
        ephem   = p.thirdBodyEphemerides.ephems[i]

        # Get gravitational parameter
        μ   = AstroUtils.getTargetGM(ephem)

        # Get position
        s   = getPosition(ephem, ts)

        # Compute relative position of spacecraft with respect to the third body
        d   = SVector(cart[1] - s[1], cart[2] - s[2], cart[3] - s[3])

        # Compute Battin's F(q) function
        Fqk = computeFk(view(cart, 1:3), s)

        # Compute acceleration
        dk      = sqrt(d[1]*d[1] + d[2]*d[2] + d[3]*d[3])
        dkInv   = 1.0 / dk
        dk3Inv  = dkInv * dkInv * dkInv
        ax     -= μ*dk3Inv*(cart[1] + Fqk*s[1])
        ay     -= μ*dk3Inv*(cart[2] + Fqk*s[2])
        az     -= μ*dk3Inv*(cart[3] + Fqk*s[3])

        # Compute partial of 1 / dk^3 w.r.t. r 
        neg3dk5Inv  = -3.0 * dk3Inv * dkInv * dkInv 
        ddk3Invdr   = SVector(neg3dk5Inv * d[1], neg3dk5Inv * d[2], neg3dk5Inv * d[3])

        # Compute partial of Fk w.r.t r
        dFkdR       = computePartialFkWrtR(view(cart, 1:3), s)

        # Compute full partial
        rPFkSk      = SVector(cart[1] + Fqk*s[1], cart[2] + Fqk*s[2], cart[3] + Fqk*s[3])
        term1       = SMatrix{3,3}(rPFkSk[1]*ddk3Invdr[1], rPFkSk[2]*ddk3Invdr[1], rPFkSk[3]*ddk3Invdr[1],
                                   rPFkSk[1]*ddk3Invdr[2], rPFkSk[2]*ddk3Invdr[2], rPFkSk[3]*ddk3Invdr[2],
                                   rPFkSk[1]*ddk3Invdr[3], rPFkSk[2]*ddk3Invdr[3], rPFkSk[3]*ddk3Invdr[3])
        term2       = SMatrix{3,3}(dk3Inv*(1.0 + s[1]*dFkdR[1]), dk3Inv*(s[2]*dFkdR[1]), dk3Inv*(s[3]*dFkdR[1]),
                                   dk3Inv*(s[1]*dFkdR[2]), dk3Inv*(1.0 + s[2]*dFkdR[2]), dk3Inv*(s[3]*dFkdR[2]),
                                   dk3Inv*(s[1]*dFkdR[3]), dk3Inv*(s[2]*dFkdR[3]), dk3Inv*(1.0 + s[3]*dFkdR[3]))

        daxdrx     -= μ*(term1[1,1] + term2[1,1])
        daxdry     -= μ*(term1[1,2] + term2[1,2])
        daxdrz     -= μ*(term1[1,3] + term2[1,3])
        daydrx     -= μ*(term1[2,1] + term2[2,1])
        daydry     -= μ*(term1[2,2] + term2[2,2])
        daydrz     -= μ*(term1[2,3] + term2[2,3])
        dazdrx     -= μ*(term1[3,1] + term2[3,1])
        dazdry     -= μ*(term1[3,2] + term2[3,2])
        dazdrz     -= μ*(term1[3,3] + term2[3,3])
    end

    return (SVector(ax, ay, az),
            SMatrix{3,7}(daxdrx, daydrx, dazdrx,
                         daxdry, daydry, dazdry,
                         daxdrz, daydrz, dazdrz,
                         daxdvx, daydvx, dazdvx,
                         daxdvy, daydvy, dazdvy,
                         daxdvz, daydvz, dazdvz,
                         daxdm,  daydm,  dazdm))
end

function computeQk(r, s)
    sds = s[1]*s[1] + s[2]*s[2] + s[3]*s[3]
    qk  = (r[1]*(r[1] - 2.0*s[1]) + 
           r[2]*(r[2] - 2.0*s[2]) + 
           r[3]*(r[3] - 2.0*s[3])) / sds
    return qk
end

function computePartialQkWrtR(r, s)
    sds = s[1]*s[1] + s[2]*s[2] + s[3]*s[3]
    twoSdsInv = 2.0 / sds 
    return SVector(twoSdsInv * (r[1] - s[1]), 
                   twoSdsInv * (r[2] - s[2]), 
                   twoSdsInv * (r[3] - s[3]))
end

function computeFk(qk)
    Fk = qk*(3.0 + 3.0*qk + qk*qk) / (1.0 + sqrt(1.0 + qk)^3)
    return Fk
end

function computeFk(r, s)
    return computeFk(computeQk(r, s))
end

function computePartialFkWrtQk(qk)
    sqrtOnePQk  = sqrt(1.0 + qk)
    return 3.0 * (qk^3*sqrtOnePQk + 2.0*(1.0 + sqrtOnePQk) + qk*qk*(2.0 + 3.0*sqrtOnePQk) + qk*(4.0 + 3.0*sqrtOnePQk)) / 
        (2.0*(1.0 + sqrtOnePQk*sqrtOnePQk*sqrtOnePQk)^2)
end

function computePartialFkWrtR(r,s)
    dFkdQk = computePartialFkWrtQk(computeQk(r, s))
    dQkdR  = computePartialQkWrtR(r, s)
    return SVector(dFkdQk*dQkdR[1], dFkdQk*dQkdR[2], dFkdQk*dQkdR[3])
end