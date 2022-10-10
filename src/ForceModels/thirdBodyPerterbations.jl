
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
        sds = s[1]*s[1] + s[2]*s[2] + s[3]*s[3]
        qk  = (cart[1]*(cart[1] - 2.0*s[1]) + 
               cart[2]*(cart[2] - 2.0*s[2]) + 
               cart[3]*(cart[3] - 2.0*s[3])) / sds
        Fqk = qk*(3.0 + 3.0*qk + qk*qk) / (1.0 + sqrt(1.0 + qk)^3) 

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