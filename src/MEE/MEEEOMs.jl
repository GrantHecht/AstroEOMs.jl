# MEE dynamics with no control out-of-place with SVector
function meeEomNoControl(u::AbstractArray, p::MEEParams, t)
    # Compute perterbing accelerations
    ap = meeComputePerterbations(u, p, t)

    # Compute requirements
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / p.mu)

    # Compute dynamics
    dp      = 2.0*u[1]*wInv*srtpbmu*ap[2]
    df      = srtpbmu*(ap[1]*sin(u[6]) + 
                       wInv*((w + 1.0)*cos(u[6]) + u[2])*ap[2] -
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3]) 
    dg      = srtpbmu*(-ap[1]*cos(u[6]) + 
                       wInv*((w + 1.0)*sin(u[6]) + u[3])*ap[2] +
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3])
    dh      = 0.5*srtpbmu*s2*wInv*cos(u[6])*ap[3]
    dk      = 0.5*srtpbmu*s2*wInv*sin(u[6])*ap[3]
    dL      = sqrt(p.mu*u[1])*wbp*wbp + 
                wInv*srtpbmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3]
    return SVector(dp,df,dg,dh,dk,dL)
end

# MEE dynamics with no control in-place
function meeEomNoControl!(du::AbstractArray, u::AbstractArray, p::MEEParams, t)
    # Compute perterbing accelerations
    ap = meeComputePerterbations(u, p, t)

    # Compute requirements
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / p.mu)

    # Compute dynamics
    du[1]   = 2.0*u[1]*wInv*srtpbmu*ap[2]
    du[2]   = srtpbmu*(ap[1]*sin(u[6]) + 
                       wInv*((w + 1.0)*cos(u[6]) + u[2])*ap[2] -
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3]) 
    du[3]   = srtpbmu*(-ap[1]*cos(u[6]) + 
                       wInv*((w + 1.0)*sin(u[6]) + u[3])*ap[2] +
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3])
    du[4]   = 0.5*srtpbmu*s2*wInv*cos(u[6])*ap[3]
    du[5]   = 0.5*srtpbmu*s2*wInv*sin(u[6])*ap[3]
    du[6]   = sqrt(p.mu*u[1])*wbp*wbp + 
                wInv*srtpbmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*ap[3]
    return nothing
end

# MEE dynamics w/ control out-of-place with SVector
function meeEomControl(u::AbstractArray, p::MEEParams, t, au::AbstractArray)
    # Compute perterbing accelerations
    ap = meeComputePerterbations(u, p, t)

    # Compute total perterbing accelerations
    apt     = SVector(ap[1] + au[1], ap[2] + au[2], ap[3] + au[3])

    # Compute requirements
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / p.mu)

    # Compute dynamics
    dp      = 2.0*u[1]*wInv*srtpbmu*apt[2]
    df      = srtpbmu*(apt[1]*sin(u[6]) + 
                       wInv*((w + 1.0)*cos(u[6]) + u[2])*apt[2] -
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3]) 
    dg      = srtpbmu*(-apt[1]*cos(u[6]) + 
                       wInv*((w + 1.0)*sin(u[6]) + u[3])*apt[2] +
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3])
    dh      = 0.5*srtpbmu*s2*wInv*cos(u[6])*apt[3]
    dk      = 0.5*srtpbmu*s2*wInv*sin(u[6])*apt[3]
    dL      = sqrt(p.mu*u[1])*wbp*wbp + 
                wInv*srtpbmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3]
    return SVector(dp,df,dg,dh,dk,dL)
end

# MEE dynamics w/ control in-place with SVector
function meeEomControl!(du::AbstractArray, u::AbstractArray, t, au::AbstractArray)
    # Compute perterbing accelerations
    ap = meeComputePerterbations(u, p, t)

    # Compute total perterbing accelerations
    apt     = SVector(ap[1] + au[1], ap[2] + au[2], ap[3] + au[3])

    # Compute requirements
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / p.mu)

    # Compute dynamics
    du[1]   = 2.0*u[1]*wInv*srtpbmu*apt[2]
    du[2]   = srtpbmu*(apt[1]*sin(u[6]) + 
                       wInv*((w + 1.0)*cos(u[6]) + u[2])*apt[2] -
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3]) 
    du[3]   = srtpbmu*(-apt[1]*cos(u[6]) + 
                       wInv*((w + 1.0)*sin(u[6]) + u[3])*apt[2] +
                       u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3])
    du[4]   = 0.5*srtpbmu*s2*wInv*cos(u[6])*apt[3]
    du[5]   = 0.5*srtpbmu*s2*wInv*sin(u[6])*apt[3]
    du[6]   = sqrt(p.mu*u[1])*wbp*wbp + 
                wInv*srtpbmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*apt[3]
    return nothing
end

function meeComputePerterbations(u::AbstractArray, p::MEEParams, t)
    Δr = 0.0
    Δt = 0.0
    Δn = 0.0
    if p.perterbations == true
        # Compute cartesian state
        cart    = convertState(u, MEE, Cartesian, p.mu)

        # Compute rotation from MEE to Inertial
        mee2I   = mee2Inertial(cart)

        # Compute gravitational perterbations
        if p.thirdBodyPerterbations == true
            # Get thirdy body perterbations in inertial frame
            atbI    = computeThirdBodyPerterbations(cart, t, p)

            # Rotate perterbations to RTN frame
            Δr     += mee2I[1,1]*atbI[1] + mee2I[2,1]*atbI[2] + mee2I[3,1]*atbI[3]
            Δt     += mee2I[1,2]*atbI[1] + mee2I[2,2]*atbI[2] + mee2I[3,2]*atbI[3]
            Δn     += mee2I[1,3]*atbI[1] + mee2I[2,3]*atbI[2] + mee2I[3,3]*atbI[3]
        end
    end
    # Return perterbing accelerations
    return SVector(Δr, Δt, Δn)
end

function mee2Inertial(cart)
    # Compute rotation from MEE to Inertial
    r     = sqrt(cart[1]*cart[1] + cart[2]*cart[2] + cart[3]*cart[3])
    hvec  = SVector(cart[2]*cart[6] - cart[3]*cart[5],
                    cart[3]*cart[4] - cart[1]*cart[6],
                    cart[1]*cart[5] - cart[2]*cart[4])
    hmag  = sqrt(hvec[1]*hvec[1] + hvec[2]*hvec[2] + hvec[3]*hvec[3])

    hrvec = SVector(hvec[2]*cart[3] - hvec[3]*cart[2],
                    hvec[3]*cart[1] - hvec[1]*cart[3],
                    hvec[1]*cart[2] - hvec[2]*cart[1])
    hrmag = sqrt(hrvec[1]*hrvec[1] + hrvec[2]*hrvec[2] + hrvec[3]*hrvec[3])

    irhat = SVector(cart[1] / r, cart[2] / r, cart[3] / r)
    inhat = SVector(hvec[1] / hmag, hvec[2] / hmag, hvec[3] / hmag)
    ithat = SVector(hrvec[1] / hrmag, hrvec[2] / hrmag, hrvec[3] / hrmag)

    return SMatrix{3,3}(irhat[1], irhat[2], irhat[3], 
                        ithat[1], ithat[2], ithat[3],
                        inhat[1], inhat[2], inhat[3])
end