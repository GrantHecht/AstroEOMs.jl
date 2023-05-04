

function computeJ2Perterbations(mee, p::AbstractEOMParams)
    # Compute requirements
    w  = 1.0 + mee[2]*cos(mee[6]) + mee[3]*sin(mee[6]) 
    s2 = 1.0 + mee[4]*mee[4] + mee[5]*mee[5]
    r  = mee[1] / w
    r4 = r*r*r*r

    # Compute acceleration due to J2
    Δr = -3.0*p.mu*p.J2*p.Rcb^2 * 
            (1.0 - 12.0*(mee[4]*sin(mee[6]) - mee[5]*cos(mee[6]))^2 / s2^2) / (2.0*r4)
    Δt = -12.0*p.mu*p.J2*p.Rcb^2 * 
            ((mee[4]*sin(mee[6]) - mee[5]*cos(mee[6])) * 
             (mee[4]*cos(mee[6]) + mee[5]*sin(mee[6])) / s2^2) / r4
    Δn = -6.0*p.mu*p.J2*p.Rcb^2 * 
            ((1.0 - mee[4]*mee[4] - mee[5]*mee[5]) * 
             (mee[4]*sin(mee[6]) - mee[5]*cos(mee[6])) / s2^2) / r4

    # Return 
    return SVector(Δr, Δt, Δn)
end