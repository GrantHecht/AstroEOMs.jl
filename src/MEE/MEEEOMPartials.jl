function meeComputeDynamicsStatePartials(u::AbstractArray, p::MEEParams, t, au)
    # Unscale inputs
    t_us    = unscaleTime(p, t)
    u_us    = unscaleStateWithMass(p, u)
    au_us   = unscaleAcceleration(p, au)

    # Compute unscaled partials and scale
    dfdx    = scaleMEEDynamicsPartials(p, 
        meeComputeDynamicsStatePartialsUnscaled(u_us, p, t_us, au_us))

    return dfdx
end

function meeComputeDynamicsTimePartials(u::AbstractArray, p::MEEParams, t)
    # Unscale inputs
    t_us    = unscaleTime(p, t)
    u_us    = unscaleStateWithMass(p, u)

    dapdt   = begin 
        # Compute unscaled acceleration partials
        dapdt_us = meeComputePerterbationTimePartials(u_us, p, t_us) 

        # Scale acceleration
        sf = p.TU*p.TU*p.TU / p.LU
        SVector(sf*dapdt_us[1], sf*dapdt_us[2], sf*dapdt_us[3])
    end

    # Compute requirements
    mu = getScaledGravityParameter(p)
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    srtpbmu = sqrt(u[1] / mu)

    # Compute B 
    B       = SMatrix{6,3}(0.0, srtpbmu*sin(u[6]), -srtpbmu*cos(u[6]), 0.0, 0.0, 0.0,
                2.0*u[1]*wInv*srtpbmu, srtpbmu*wInv*((w + 1.0)*cos(u[6]) + u[2]), srtpbmu*wInv*((w + 1.0)*sin(u[6]) + u[3]), 0.0, 0.0, 0.0,
                0.0, -srtpbmu*u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), srtpbmu*u[2]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), 
                0.5*srtpbmu*s2*wInv*cos(u[6]), 0.5*srtpbmu*s2*wInv*sin(u[6]), srtpbmu*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])))

    return SVector(B[1,2]*dapdt[2],
                   B[2,1]*dapdt[1] + B[2,2]*dapdt[2] + B[2,3]*dapdt[3],
                   B[3,1]*dapdt[1] + B[3,2]*dapdt[2] + B[3,3]*dapdt[3],
                   B[4,3]*dapdt[3],
                   B[5,3]*dapdt[3],
                   B[6,3]*dapdt[3],
                   0.0)
end

function meeComputeDynamicsThrustScaleFactorPartials(u::AbstractArray, p::Tuple{SimpleSpacecraft,MEEParams}, t, α)
    # Compute partial of acceleration due to thrust w.r.t. scale factor
    daumagdsf = p[1].tMax*p[2].TU*p[2].TU / (1000.0 * p[2].LU * p[2].MU * u[7])
    daudsf    = SVector(daumagdsf*α[1], daumagdsf*α[2], daumagdsf*α[3]) 

    # Compute requirements
    mu      = getScaledGravityParameter(p[2])
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / mu)
    c       = p[1].c * p[2].TU / (1000.0 * p[2].LU)
    cInv    = 1.0 / c

    # Compute B 
    B       = SMatrix{7,3}(0.0, srtpbmu*sin(u[6]), -srtpbmu*cos(u[6]), 0.0, 0.0, 0.0, -u[7]*cInv,
                2.0*u[1]*wInv*srtpbmu, srtpbmu*wInv*((w + 1.0)*cos(u[6]) + u[2]), srtpbmu*wInv*((w + 1.0)*sin(u[6]) + u[3]), 0.0, 0.0, 0.0, -u[7]*cInv,
                0.0, -srtpbmu*u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), srtpbmu*u[2]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), 
                0.5*srtpbmu*s2*wInv*cos(u[6]), 0.5*srtpbmu*s2*wInv*sin(u[6]), srtpbmu*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), -u[7]*cInv)

    # Compute and return partial
    return SVector(B[1,1]*daudsf[1] + B[1,2]*daudsf[2] + B[1,3]*daudsf[3],
                   B[2,1]*daudsf[1] + B[2,2]*daudsf[2] + B[2,3]*daudsf[3],
                   B[3,1]*daudsf[1] + B[3,2]*daudsf[2] + B[3,3]*daudsf[3],
                   B[4,1]*daudsf[1] + B[4,2]*daudsf[2] + B[4,3]*daudsf[3],
                   B[5,1]*daudsf[1] + B[5,2]*daudsf[2] + B[5,3]*daudsf[3],
                   B[6,1]*daudsf[1] + B[6,2]*daudsf[2] + B[6,3]*daudsf[3],
                   B[7,1]*daudsf[1] + B[7,2]*daudsf[2] + B[7,3]*daudsf[3])
end

function meeComputeDynamicsThrustUnitVectorPartials(u::AbstractArray, p::Tuple{SimpleSpacecraft,MEEParams}, t, sf)
    # Compute partial of acceleration due to thrust w.r.t. scale factor
    auMag   = sf*p[1].tMax*p[2].TU*p[2].TU / (1000.0 * p[2].LU * p[2].MU * u[7])
    daudα   = SVector(auMag, auMag, auMag) 

    # Compute requirements
    mu      = getScaledGravityParameter(p[2])
    s2      = 1.0 + u[4]*u[4] + u[5]*u[5]
    w       = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    wInv    = 1.0 / w
    wbp     = w / u[1]
    srtpbmu = sqrt(u[1] / mu)
    c       = p[1].c * p[2].TU / (1000.0 * p[2].LU)
    cInv    = 1.0 / c

    # Compute B 
    B       = SMatrix{7,3}(0.0, srtpbmu*sin(u[6]), -srtpbmu*cos(u[6]), 0.0, 0.0, 0.0, -u[7]*cInv,
                2.0*u[1]*wInv*srtpbmu, srtpbmu*wInv*((w + 1.0)*cos(u[6]) + u[2]), srtpbmu*wInv*((w + 1.0)*sin(u[6]) + u[3]), 0.0, 0.0, 0.0, -u[7]*cInv,
                0.0, -srtpbmu*u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), srtpbmu*u[2]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), 
                0.5*srtpbmu*s2*wInv*cos(u[6]), 0.5*srtpbmu*s2*wInv*sin(u[6]), srtpbmu*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), -u[7]*cInv)

    # Compute and return partial
    return SMatrix{7,3}(B[1,1]*daudα[1], B[2,1]*daudα[1], B[3,1]*daudα[1], B[4,1]*daudα[1], B[5,1]*daudα[1], B[6,1]*daudα[1], B[7,1]*daudα[1],
                        B[1,2]*daudα[2], B[2,2]*daudα[2], B[3,2]*daudα[1], B[4,2]*daudα[2], B[5,2]*daudα[2], B[6,2]*daudα[2], B[7,2]*daudα[2],
                        B[1,3]*daudα[3], B[2,3]*daudα[3], B[3,3]*daudα[1], B[4,3]*daudα[3], B[5,3]*daudα[3], B[6,3]*daudα[3], B[7,3]*daudα[3])
end

# All inputs are unscaled
function meeComputeDynamicsStatePartialsUnscaled(u::AbstractArray, p::MEEParams, t, au)
    Δ, dΔdmee = begin # Get total perterbations and partials
        # Compute perterbations and partials
        ap, dapdmee = meeComputePerterbationStatePartials(u, p, t) 

        # Compute total perterbation (i.e. celestial phenomena + thrust)
        Δ = SVector(ap[1] + au[1], ap[2] + au[2], ap[3] + au[3])

        # Compute total perterbation partials
        mInv   = 1.0 / u[7]
        dΔdmee = SMatrix{3,7}(dapdmee[1,1], dapdmee[2,1], dapdmee[3,1],
                              dapdmee[1,2], dapdmee[2,2], dapdmee[3,2],
                              dapdmee[1,3], dapdmee[2,3], dapdmee[3,3],
                              dapdmee[1,4], dapdmee[2,4], dapdmee[3,4],
                              dapdmee[1,5], dapdmee[2,5], dapdmee[3,5],
                              dapdmee[1,6], dapdmee[2,6], dapdmee[3,6],
                              dapdmee[1,7] - mInv*au[1], 
                              dapdmee[2,7] - mInv*au[2], 
                              dapdmee[3,7] - mInv*au[3]) 
        (Δ,dΔdmee)
    end

    # Compute requirements
    ss = 1.0 + u[4]*u[4] + u[5]*u[5]
    w  = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])

    wInv    = 1.0 / w
    pInv    = 1.0 / u[1]
    sqrtPmu = sqrt(u[1] / p.mu)

    # Partials of pdot 
    fpt     = 2.0*u[1]*wInv*sqrtPmu
    dfpdp   = 3.0*wInv*sqrtPmu*Δ[2] + fpt*dΔdmee[2,1]
    dfpdf   = -2.0*u[1]*wInv*wInv*sqrtPmu*cos(u[6])*Δ[2] + fpt*dΔdmee[2,2]
    dfpdg   = -2.0*u[1]*wInv*wInv*sqrtPmu*sin(u[6])*Δ[2] + fpt*dΔdmee[2,3]
    dfpdh   = fpt*dΔdmee[2,4]
    dfpdk   = fpt*dΔdmee[2,5]
    dfpdL   = -2.0*u[1]*wInv*wInv*sqrtPmu*(u[3]*cos(u[6]) - u[2]*sin(u[6]))*Δ[2] + fpt*dΔdmee[2,6]
    dfpdm   = fpt*dΔdmee[2,7]

    # Partials of fdot
    ffr     =  sqrtPmu*sin(u[6])
    fft     =  sqrtPmu*wInv*((w + 1.0)*cos(u[6]) + u[2])
    ffn     = -sqrtPmu*wInv*u[3]*(u[4]*sin(u[6]) - u[5]*cos(u[6]))
    dffdp   = 0.5*pInv*wInv*sqrtPmu*(w*sin(u[6])*Δ[1] + u[3]*(u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3] + 
                (u[2] + u[2]*cos(u[6])^2 + cos(u[6])*(2.0 + u[3]*sin(u[6])))*Δ[2]) + 
                ffr*dΔdmee[1,1] + fft*dΔdmee[2,1] + ffn*dΔdmee[3,1]
    dffdf   = wInv*wInv*sqrtPmu*(sin(u[6])*(u[3] + sin(u[6]))*Δ[2] - u[3]*cos(u[6])*(u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3]) + 
                ffr*dΔdmee[1,2] + fft*dΔdmee[2,2] + ffn*dΔdmee[3,2]
    dffdg   = wInv*wInv*sqrtPmu*(-(u[2] + cos(u[6]))*sin(u[6])*Δ[2] + (1.0 + u[2]*cos(u[6]))*(u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3]) + 
                ffr*dΔdmee[1,3] + fft*dΔdmee[2,3] + ffn*dΔdmee[3,3]
    dffdh   = -u[3]*wInv*sqrtPmu*sin(u[6])*Δ[3] + 
                ffr*dΔdmee[1,4] + fft*dΔdmee[2,4] + ffn*dΔdmee[3,4]
    dffdk   = u[3]*wInv*sqrtPmu*cos(u[6])*Δ[3] + 
                ffr*dΔdmee[1,5] + fft*dΔdmee[2,5] + ffn*dΔdmee[3,5]
    dffdL   = -wInv*wInv*sqrtPmu*(-w*w*cos(u[6])*Δ[1] + u[3]*(u[2]*u[4] + u[3]*u[5] + u[4]*cos(u[6]) + u[5]*sin(u[6]))*Δ[3] + 
                (u[2]*u[3]*cos(u[6]) + 3.0*u[3]*sin(u[6])^2 + u[3]*u[3]*sin(u[6])^3 + cos(u[6])^2*(u[3] + u[2]*u[2]*sin(u[6])) + 
                 u[2]*sin(2.0*u[6]) + sin(u[6])*(2.0 - u[2]*u[2] + u[2]*u[3]*sin(2.0*u[6])))*Δ[2]) + 
                ffr*dΔdmee[1,6] + fft*dΔdmee[2,6] + ffn*dΔdmee[3,6]
    dffdm   = ffr*dΔdmee[1,7] + fft*dΔdmee[2,7] + ffn*dΔdmee[3,7]

    # Partials of gdot
    fgr     = -sqrtPmu*cos(u[6])
    fgt     =  sqrtPmu*wInv*((w + 1.0)*sin(u[6]) + u[3])
    fgn     =  sqrtPmu*wInv*u[2]*(u[4]*sin(u[6]) - u[5]*cos(u[6]))
    dfgdp   = 0.5*pInv*sqrtPmu*(-cos(u[6])*Δ[1] + wInv*(u[2]*(-u[5]*cos(u[6]) + u[4]*sin(u[6]))*Δ[3] + 
                (u[3] + (2.0 + u[2]*cos(u[6]))*sin(u[6]) + u[3]*sin(u[6])^2)*Δ[2])) + 
                fgr*dΔdmee[1,1] + fgt*dΔdmee[2,1] + fgn*dΔdmee[3,1]
    dfgdf   = wInv*wInv*sqrtPmu*(Δ[3]*u[4]*sin(u[6])*(1.0 + u[3]*sin(u[6])) - cos(u[6])*(Δ[2]*u[3] + Δ[3]*u[5] + (Δ[2] + Δ[3]*u[3]*u[5])*sin(u[6]))) +
                fgr*dΔdmee[1,2] + fgt*dΔdmee[2,2] + fgn*dΔdmee[3,2]
    dfgdg   = wInv*wInv*sqrtPmu*(Δ[2] - (Δ[2] + Δ[3]*u[2]*u[4])*sin(u[6])^2 + u[2]*cos(u[6])*(Δ[2] + Δ[3]*u[5]*sin(u[6]))) +
                fgr*dΔdmee[1,3] + fgt*dΔdmee[2,3] + fgn*dΔdmee[3,3]
    dfgdh   = u[2]*wInv*sqrtPmu*sin(u[6])*Δ[3] + 
                fgr*dΔdmee[1,4] + fgt*dΔdmee[2,4] + fgn*dΔdmee[3,4]
    dfgdk   = -u[2]*wInv*sqrtPmu*cos(u[6])*Δ[3] + 
                fgr*dΔdmee[1,5] + fgt*dΔdmee[2,5] + fgn*dΔdmee[3,5]
    dfgdL   = wInv*wInv*sqrtPmu*(w*w*sin(u[6])*Δ[1] + u[2]*(u[2]*u[4] + u[3]*u[5] + u[4]*cos(u[6]) + u[5]*sin(u[6]))*Δ[3] +
                (u[2]*u[2]*cos(u[6])^3 + u[2]*sin(u[6])*(u[3] + sin(u[6])) + u[2]*cos(u[6])^2*(3.0 + 2.0*u[3]*sin(u[6])) + 
                cos(u[6])*(2.0 - u[3]*u[3] + 2.0*u[3]*sin(u[6]) + u[3]*u[3]*sin(u[6])^2))*Δ[2]) + 
                fgr*dΔdmee[1,6] + fgt*dΔdmee[2,6] + fgn*dΔdmee[3,6]
    dfgdm   = fgr*dΔdmee[1,7] + fgt*dΔdmee[2,7] + fgn*dΔdmee[3,7]

    # Partials of hdot
    fhn     =  0.5*sqrtPmu*ss*wInv*cos(u[6])
    dfhdp   =  0.25*ss*pInv*wInv*sqrtPmu*cos(u[6])*Δ[3] + fhn*dΔdmee[3,1]
    dfhdf   = -0.5*ss*wInv*wInv*sqrtPmu*cos(u[6])^2*Δ[3] + fhn*dΔdmee[3,2]
    dfhdg   = -0.25*ss*wInv*wInv*sqrtPmu*sin(2.0*u[6])*Δ[3] + fhn*dΔdmee[3,3]
    dfhdh   =  u[4]*wInv*sqrtPmu*cos(u[6])*Δ[3] + fhn*dΔdmee[3,4]
    dfhdk   =  u[5]*wInv*sqrtPmu*cos(u[6])*Δ[3] + fhn*dΔdmee[3,5]
    dfhdL   = -0.5*ss*wInv*wInv*sqrtPmu*(u[3] + sin(u[6]))*Δ[3] + fhn*dΔdmee[3,6]
    dfhdm   =  fhn*dΔdmee[3,7]

    # Partials of kdot
    fkn     =  0.5*sqrtPmu*ss*wInv*sin(u[6])
    dfkdp   =  0.25*ss*pInv*wInv*sqrtPmu*sin(u[6])*Δ[3] + fkn*dΔdmee[3,1]
    dfkdf   = -0.25*ss*wInv*wInv*sqrtPmu*sin(2.0*u[6])*Δ[3] + fkn*dΔdmee[3,2]
    dfkdg   = -0.5*ss*wInv*wInv*sqrtPmu*sin(u[6])^2*Δ[3] + fkn*dΔdmee[3,3]
    dfkdh   =  u[4]*wInv*sqrtPmu*sin(u[6])*Δ[3] + fkn*dΔdmee[3,4]
    dfkdk   =  u[5]*wInv*sqrtPmu*sin(u[6])*Δ[3] + fkn*dΔdmee[3,5]
    dfkdL   =  0.5*ss*wInv*wInv*sqrtPmu*(u[2] + cos(u[6]))*Δ[3] + fkn*dΔdmee[3,6]
    dfkdm   =  fkn*dΔdmee[3,7]

    # Partials of Ldot
    fLn     = wInv*sqrtPmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))
    dfLdp   = 0.5*wInv*(-3.0*pInv*pInv*pInv*sqrt(p.mu*u[1])*w*w*w + pInv*sqrtPmu*(u[4]*sin(u[6]) - u[5]*cos(u[6]))*Δ[3]) +
                fLn*dΔdmee[3,1]
    dfLdf   = 2.0*w*pInv*pInv*sqrt(p.mu*u[1])*cos(u[6]) + wInv*wInv*sqrtPmu*cos(u[6])*(u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3] +
                fLn*dΔdmee[3,2]
    dfLdg   = 2.0*w*pInv*pInv*sqrt(p.mu*u[1])*sin(u[6]) + wInv*wInv*sqrtPmu*sin(u[6])*(u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3] + 
                fLn*dΔdmee[3,3]
    dfLdh   = wInv*sqrtPmu*sin(u[6])*Δ[3] + fLn*dΔdmee[3,4]
    dfLdk   = -wInv*sqrtPmu*cos(u[6])*Δ[3] + fLn*dΔdmee[3,5]
    dfLdL   = 2.0*w*pInv*pInv*sqrt(p.mu*u[1])*(u[3]*cos(u[6]) - u[2]*sin(u[6])) + wInv*wInv*sqrtPmu*(u[3]*cos(u[6]) - u[2]*sin(u[6])) * 
                (u[5]*cos(u[6]) - u[4]*sin(u[6]))*Δ[3] + wInv*sqrtPmu*(u[4]*cos(u[6]) + u[5]*sin(u[6]))*Δ[3] +
                fLn*dΔdmee[3,6]
    dfLdm   = fLn*dΔdmee[3,7]

    # Partials of mdot
    dfmdp   = 0.0
    dfmdf   = 0.0
    dfmdg   = 0.0
    dfmdh   = 0.0
    dfmdk   = 0.0
    dfmdL   = 0.0
    dfmdm   = 0.0

    return SMatrix{7,7}(dfpdp, dffdp, dfgdp, dfhdp, dfkdp, dfLdp, dfmdp,
                        dfpdf, dffdf, dfgdf, dfhdf, dfkdf, dfLdf, dfmdf,
                        dfpdg, dffdg, dfgdg, dfhdg, dfkdg, dfLdg, dfmdg,
                        dfpdh, dffdh, dfgdh, dfhdh, dfkdh, dfLdh, dfmdh,
                        dfpdk, dffdk, dfgdk, dfhdk, dfkdk, dfLdk, dfmdk,
                        dfpdL, dffdL, dfgdL, dfhdL, dfkdL, dfLdL, dfmdL,
                        dfpdm, dffdm, dfgdm, dfhdm, dfkdm, dfLdm, dfmdm)
end

function meeComputePerterbationTimePartials(u::AbstractArray, p::MEEParams, t)
    Δrdt = 0.0
    Δtdt = 0.0
    Δndt = 0.0
    if p.perterbations == true
        # Compute cartesian state and partial of cartesian state w.r.t. mee state
        cart        = convertState(view(u, 1:6), MEE, Cartesian, p.mu)

        # Construct cartesian state with mass
        cartWithMass    = SVector(cart[1], cart[2], cart[3], cart[4],
                            cart[5], cart[6], u[7])

        # Compute partials of rotation matrix w.r.t cartesian state
        mee2I   = mee2Inertial(cart)

        if p.thirdBodyPerterbations == true
            # Compute 3rd body perterbation partials
            datbidt = computeThirdBodyPerterbationTimePartials(cartWithMass, t, p)

            # Rotate perterbations to RTN frame
            Δrdt   += mee2I[1,1]*datbidt[1] + mee2I[2,1]*datbidt[2] + mee2I[3,1]*datbidt[3]
            Δtdt   += mee2I[1,2]*datbidt[1] + mee2I[2,2]*datbidt[2] + mee2I[3,2]*datbidt[3]
            Δndt   += mee2I[1,3]*datbidt[1] + mee2I[2,3]*datbidt[2] + mee2I[3,3]*datbidt[3]
        end
    end
    return SVector(Δrdt, Δtdt, Δndt)
end

function meeComputePerterbationStatePartials(u::AbstractArray, p::MEEParams, t)
    Δr      = 0.0
    Δt      = 0.0
    Δn      = 0.0
    Δrdp    = 0.0
    Δrdf    = 0.0
    Δrdg    = 0.0
    Δrdh    = 0.0
    Δrdk    = 0.0
    ΔrdL    = 0.0
    Δrdm    = 0.0
    Δtdp    = 0.0
    Δtdf    = 0.0
    Δtdg    = 0.0
    Δtdh    = 0.0
    Δtdk    = 0.0
    ΔtdL    = 0.0
    Δtdm    = 0.0
    Δndp    = 0.0
    Δndf    = 0.0
    Δndg    = 0.0
    Δndh    = 0.0
    Δndk    = 0.0
    ΔndL    = 0.0
    Δndm    = 0.0
    if p.perterbations == true
        # Compute cartesian state and partial of cartesian state w.r.t. mee state
        cart        = convertState(view(u, 1:6), MEE, Cartesian, p.mu)
        dcartdmee   = convertStatePartials(u, MEE, Cartesian, p.mu)

        # Construct cartesian state with mass
        cartWithMass    = SVector(cart[1], cart[2], cart[3], cart[4],
                            cart[5], cart[6], u[7])

        # Compute partials of rotation matrix w.r.t cartesian state
        (s1dAdcart, s2dAdcart, s3dAdcart) = mee2InertialPartialsWrtCartesianState(cart)
        mee2I   = mee2Inertial(cart)

        if p.thirdBodyPerterbations == true
            # Compute 3rd body perterbation partials
            atbi, datbidcart = computeThirdBodyPerterbationsAndPartials(cartWithMass, t, p)

            # Rotate perterbations to RTN frame
            Δr     += mee2I[1,1]*atbi[1] + mee2I[2,1]*atbi[2] + mee2I[3,1]*atbi[3]
            Δt     += mee2I[1,2]*atbi[1] + mee2I[2,2]*atbi[2] + mee2I[3,2]*atbi[3]
            Δn     += mee2I[1,3]*atbi[1] + mee2I[2,3]*atbi[2] + mee2I[3,3]*atbi[3]

            # Compute partial of atbRTN w.r.t cartesian state (no mass partial for term1 because corresponding terms are zero)
            term1   = SMatrix{3,6}(atbi[1]*s1dAdcart[1,1] + atbi[2]*s1dAdcart[2,1] + atbi[3]*s1dAdcart[3,1],
                                   atbi[1]*s2dAdcart[1,1] + atbi[2]*s2dAdcart[2,1] + atbi[3]*s2dAdcart[3,1],
                                   atbi[1]*s3dAdcart[1,1] + atbi[2]*s3dAdcart[2,1] + atbi[3]*s3dAdcart[3,1],
                                   atbi[1]*s1dAdcart[1,2] + atbi[2]*s1dAdcart[2,2] + atbi[3]*s1dAdcart[3,2],
                                   atbi[1]*s2dAdcart[1,2] + atbi[2]*s2dAdcart[2,2] + atbi[3]*s2dAdcart[3,2],
                                   atbi[1]*s3dAdcart[1,2] + atbi[2]*s3dAdcart[2,2] + atbi[3]*s3dAdcart[3,2],
                                   atbi[1]*s1dAdcart[1,3] + atbi[2]*s1dAdcart[2,3] + atbi[3]*s1dAdcart[3,3],
                                   atbi[1]*s2dAdcart[1,3] + atbi[2]*s2dAdcart[2,3] + atbi[3]*s2dAdcart[3,3],
                                   atbi[1]*s3dAdcart[1,3] + atbi[2]*s3dAdcart[2,3] + atbi[3]*s3dAdcart[3,3],
                                   atbi[1]*s1dAdcart[1,4] + atbi[2]*s1dAdcart[2,4] + atbi[3]*s1dAdcart[3,4],
                                   atbi[1]*s2dAdcart[1,4] + atbi[2]*s2dAdcart[2,4] + atbi[3]*s2dAdcart[3,4],
                                   atbi[1]*s3dAdcart[1,4] + atbi[2]*s3dAdcart[2,4] + atbi[3]*s3dAdcart[3,4],
                                   atbi[1]*s1dAdcart[1,5] + atbi[2]*s1dAdcart[2,5] + atbi[3]*s1dAdcart[3,5],
                                   atbi[1]*s2dAdcart[1,5] + atbi[2]*s2dAdcart[2,5] + atbi[3]*s2dAdcart[3,5],
                                   atbi[1]*s3dAdcart[1,5] + atbi[2]*s3dAdcart[2,5] + atbi[3]*s3dAdcart[3,5],
                                   atbi[1]*s1dAdcart[1,6] + atbi[2]*s1dAdcart[2,6] + atbi[3]*s1dAdcart[3,6],
                                   atbi[1]*s2dAdcart[1,6] + atbi[2]*s2dAdcart[2,6] + atbi[3]*s2dAdcart[3,6],
                                   atbi[1]*s3dAdcart[1,6] + atbi[2]*s3dAdcart[2,6] + atbi[3]*s3dAdcart[3,6])

            term2   = SMatrix{3,7}(mee2I[1,1]*datbidcart[1,1] + mee2I[2,1]*datbidcart[2,1] + mee2I[3,1]*datbidcart[3,1],
                                   mee2I[1,2]*datbidcart[1,1] + mee2I[2,2]*datbidcart[2,1] + mee2I[3,2]*datbidcart[3,1],
                                   mee2I[1,3]*datbidcart[1,1] + mee2I[2,3]*datbidcart[2,1] + mee2I[3,3]*datbidcart[3,1],
                                   mee2I[1,1]*datbidcart[1,2] + mee2I[2,1]*datbidcart[2,2] + mee2I[3,1]*datbidcart[3,2],
                                   mee2I[1,2]*datbidcart[1,2] + mee2I[2,2]*datbidcart[2,2] + mee2I[3,2]*datbidcart[3,2],
                                   mee2I[1,3]*datbidcart[1,2] + mee2I[2,3]*datbidcart[2,2] + mee2I[3,3]*datbidcart[3,2],
                                   mee2I[1,1]*datbidcart[1,3] + mee2I[2,1]*datbidcart[2,3] + mee2I[3,1]*datbidcart[3,3],
                                   mee2I[1,2]*datbidcart[1,3] + mee2I[2,2]*datbidcart[2,3] + mee2I[3,2]*datbidcart[3,3],
                                   mee2I[1,3]*datbidcart[1,3] + mee2I[2,3]*datbidcart[2,3] + mee2I[3,3]*datbidcart[3,3],
                                   mee2I[1,1]*datbidcart[1,4] + mee2I[2,1]*datbidcart[2,4] + mee2I[3,1]*datbidcart[3,4],
                                   mee2I[1,2]*datbidcart[1,4] + mee2I[2,2]*datbidcart[2,4] + mee2I[3,2]*datbidcart[3,4],
                                   mee2I[1,3]*datbidcart[1,4] + mee2I[2,3]*datbidcart[2,4] + mee2I[3,3]*datbidcart[3,4],
                                   mee2I[1,1]*datbidcart[1,5] + mee2I[2,1]*datbidcart[2,5] + mee2I[3,1]*datbidcart[3,5],
                                   mee2I[1,2]*datbidcart[1,5] + mee2I[2,2]*datbidcart[2,5] + mee2I[3,2]*datbidcart[3,5],
                                   mee2I[1,3]*datbidcart[1,5] + mee2I[2,3]*datbidcart[2,5] + mee2I[3,3]*datbidcart[3,5],
                                   mee2I[1,1]*datbidcart[1,6] + mee2I[2,1]*datbidcart[2,6] + mee2I[3,1]*datbidcart[3,6],
                                   mee2I[1,2]*datbidcart[1,6] + mee2I[2,2]*datbidcart[2,6] + mee2I[3,2]*datbidcart[3,6],
                                   mee2I[1,3]*datbidcart[1,6] + mee2I[2,3]*datbidcart[2,6] + mee2I[3,3]*datbidcart[3,6],
                                   mee2I[1,1]*datbidcart[1,7] + mee2I[2,1]*datbidcart[2,7] + mee2I[3,1]*datbidcart[3,7],
                                   mee2I[1,2]*datbidcart[1,7] + mee2I[2,2]*datbidcart[2,7] + mee2I[3,2]*datbidcart[3,7],
                                   mee2I[1,3]*datbidcart[1,7] + mee2I[2,3]*datbidcart[2,7] + mee2I[3,3]*datbidcart[3,7])

            dartndcart  = SMatrix{3,7}(term1[1,1] + term2[1,1], term1[2,1] + term2[2,1],  term1[3,1] + term2[3,1],
                                       term1[1,2] + term2[1,2], term1[2,2] + term2[2,2],  term1[3,2] + term2[3,2],
                                       term1[1,3] + term2[1,3], term1[2,3] + term2[2,3],  term1[3,3] + term2[3,3],
                                       term1[1,4] + term2[1,4], term1[2,4] + term2[2,4],  term1[3,4] + term2[3,4],
                                       term1[1,5] + term2[1,5], term1[2,5] + term2[2,5],  term1[3,5] + term2[3,5],
                                       term1[1,6] + term2[1,6], term1[2,6] + term2[2,6],  term1[3,6] + term2[3,6],
                                       term2[1,7],              term2[2,7],               term2[3,7])

            Δrdp   += dartndcart[1,1]*dcartdmee[1,1] + dartndcart[1,2]*dcartdmee[2,1] + dartndcart[1,3]*dcartdmee[3,1] + dartndcart[1,4]*dcartdmee[4,1] + dartndcart[1,5]*dcartdmee[5,1] + dartndcart[1,6]*dcartdmee[6,1]
            Δtdp   += dartndcart[2,1]*dcartdmee[1,1] + dartndcart[2,2]*dcartdmee[2,1] + dartndcart[2,3]*dcartdmee[3,1] + dartndcart[2,4]*dcartdmee[4,1] + dartndcart[2,5]*dcartdmee[5,1] + dartndcart[2,6]*dcartdmee[6,1]
            Δndp   += dartndcart[3,1]*dcartdmee[1,1] + dartndcart[3,2]*dcartdmee[2,1] + dartndcart[3,3]*dcartdmee[3,1] + dartndcart[3,4]*dcartdmee[4,1] + dartndcart[3,5]*dcartdmee[5,1] + dartndcart[3,6]*dcartdmee[6,1]
            Δrdf   += dartndcart[1,1]*dcartdmee[1,2] + dartndcart[1,2]*dcartdmee[2,2] + dartndcart[1,3]*dcartdmee[3,2] + dartndcart[1,4]*dcartdmee[4,2] + dartndcart[1,5]*dcartdmee[5,2] + dartndcart[1,6]*dcartdmee[6,2]
            Δtdf   += dartndcart[2,1]*dcartdmee[1,2] + dartndcart[2,2]*dcartdmee[2,2] + dartndcart[2,3]*dcartdmee[3,2] + dartndcart[2,4]*dcartdmee[4,2] + dartndcart[2,5]*dcartdmee[5,2] + dartndcart[2,6]*dcartdmee[6,2]
            Δndf   += dartndcart[3,1]*dcartdmee[1,2] + dartndcart[3,2]*dcartdmee[2,2] + dartndcart[3,3]*dcartdmee[3,2] + dartndcart[3,4]*dcartdmee[4,2] + dartndcart[3,5]*dcartdmee[5,2] + dartndcart[3,6]*dcartdmee[6,2]
            Δrdg   += dartndcart[1,1]*dcartdmee[1,3] + dartndcart[1,2]*dcartdmee[2,3] + dartndcart[1,3]*dcartdmee[3,3] + dartndcart[1,4]*dcartdmee[4,3] + dartndcart[1,5]*dcartdmee[5,3] + dartndcart[1,6]*dcartdmee[6,3]
            Δtdg   += dartndcart[2,1]*dcartdmee[1,3] + dartndcart[2,2]*dcartdmee[2,3] + dartndcart[2,3]*dcartdmee[3,3] + dartndcart[2,4]*dcartdmee[4,3] + dartndcart[2,5]*dcartdmee[5,3] + dartndcart[2,6]*dcartdmee[6,3]
            Δndg   += dartndcart[3,1]*dcartdmee[1,3] + dartndcart[3,2]*dcartdmee[2,3] + dartndcart[3,3]*dcartdmee[3,3] + dartndcart[3,4]*dcartdmee[4,3] + dartndcart[3,5]*dcartdmee[5,3] + dartndcart[3,6]*dcartdmee[6,3]
            Δrdh   += dartndcart[1,1]*dcartdmee[1,4] + dartndcart[1,2]*dcartdmee[2,4] + dartndcart[1,3]*dcartdmee[3,4] + dartndcart[1,4]*dcartdmee[4,4] + dartndcart[1,5]*dcartdmee[5,4] + dartndcart[1,6]*dcartdmee[6,4]
            Δtdh   += dartndcart[2,1]*dcartdmee[1,4] + dartndcart[2,2]*dcartdmee[2,4] + dartndcart[2,3]*dcartdmee[3,4] + dartndcart[2,4]*dcartdmee[4,4] + dartndcart[2,5]*dcartdmee[5,4] + dartndcart[2,6]*dcartdmee[6,4]
            Δndh   += dartndcart[3,1]*dcartdmee[1,4] + dartndcart[3,2]*dcartdmee[2,4] + dartndcart[3,3]*dcartdmee[3,4] + dartndcart[3,4]*dcartdmee[4,4] + dartndcart[3,5]*dcartdmee[5,4] + dartndcart[3,6]*dcartdmee[6,4]
            Δrdk   += dartndcart[1,1]*dcartdmee[1,5] + dartndcart[1,2]*dcartdmee[2,5] + dartndcart[1,3]*dcartdmee[3,5] + dartndcart[1,4]*dcartdmee[4,5] + dartndcart[1,5]*dcartdmee[5,5] + dartndcart[1,6]*dcartdmee[6,5]
            Δtdk   += dartndcart[2,1]*dcartdmee[1,5] + dartndcart[2,2]*dcartdmee[2,5] + dartndcart[2,3]*dcartdmee[3,5] + dartndcart[2,4]*dcartdmee[4,5] + dartndcart[2,5]*dcartdmee[5,5] + dartndcart[2,6]*dcartdmee[6,5]
            Δndk   += dartndcart[3,1]*dcartdmee[1,5] + dartndcart[3,2]*dcartdmee[2,5] + dartndcart[3,3]*dcartdmee[3,5] + dartndcart[3,4]*dcartdmee[4,5] + dartndcart[3,5]*dcartdmee[5,5] + dartndcart[3,6]*dcartdmee[6,5]
            ΔrdL   += dartndcart[1,1]*dcartdmee[1,6] + dartndcart[1,2]*dcartdmee[2,6] + dartndcart[1,3]*dcartdmee[3,6] + dartndcart[1,4]*dcartdmee[4,6] + dartndcart[1,5]*dcartdmee[5,6] + dartndcart[1,6]*dcartdmee[6,6]
            ΔtdL   += dartndcart[2,1]*dcartdmee[1,6] + dartndcart[2,2]*dcartdmee[2,6] + dartndcart[2,3]*dcartdmee[3,6] + dartndcart[2,4]*dcartdmee[4,6] + dartndcart[2,5]*dcartdmee[5,6] + dartndcart[2,6]*dcartdmee[6,6]
            ΔndL   += dartndcart[3,1]*dcartdmee[1,6] + dartndcart[3,2]*dcartdmee[2,6] + dartndcart[3,3]*dcartdmee[3,6] + dartndcart[3,4]*dcartdmee[4,6] + dartndcart[3,5]*dcartdmee[5,6] + dartndcart[3,6]*dcartdmee[6,6]
            Δrdm   += dartndcart[1,7]
            Δtdm   += dartndcart[2,7]
            Δndm   += dartndcart[3,7]
        end
    end
    return (SVector(Δr, Δt, Δn),
            SMatrix{3,7}(Δrdp, Δtdp, Δndp,
                         Δrdf, Δtdf, Δndf,
                         Δrdg, Δtdg, Δndg,
                         Δrdh, Δtdh, Δndh,
                         Δrdk, Δtdk, Δndk,
                         ΔrdL, ΔtdL, ΔndL,
                         Δrdm, Δtdm, Δndm))
end

function mee2InertialPartialsWrtCartesianState(cart)
    # Compute requirements
    r       = sqrt(cart[1]*cart[1] + cart[2]*cart[2] + cart[3]*cart[3])
    invR    = 1.0 / r
    invR3   = invR*invR*invR

    # Compute slice 1 (∂A[:,1] / ∂cart or ∂irhat / ∂cart)
    s1  = SMatrix{3,6}((cart[2]*cart[2] + cart[3]*cart[3]) * invR3,
                       -(cart[1]*cart[2]) * invR3,
                       -(cart[1]*cart[3]) * invR3,
                       -(cart[1]*cart[2]) * invR3,
                       (cart[1]*cart[1] + cart[3]*cart[3]) * invR3,
                       -(cart[2]*cart[3]) * invR3,
                       -(cart[1]*cart[3]) * invR3,
                       -(cart[2]*cart[3]) * invR3,
                       (cart[1]*cart[1] + cart[2]*cart[2]) * invR3,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Compute slice 3 (∂A[:,3] / ∂cart or ∂inhat / ∂cart)
    s3 = begin
        hvec        = SVector(cart[2]*cart[6] - cart[3]*cart[5],
                              cart[3]*cart[4] - cart[1]*cart[6],
                              cart[1]*cart[5] - cart[2]*cart[4])
        hmag        = sqrt(hvec[1]*hvec[1] + hvec[2]*hvec[2] + hvec[3]*hvec[3])
        invHmag     = 1.0 / hmag
        invHmag3    = invHmag * invHmag * invHmag

        dhhatdh     = SMatrix{3,3}((hvec[2]*hvec[2] + hvec[3]*hvec[3]) * invHmag3,
                                -(hvec[1]*hvec[2]) * invHmag3,
                                -(hvec[1]*hvec[3]) * invHmag3,
                                -(hvec[1]*hvec[2]) * invHmag3,
                                (hvec[1]*hvec[1] + hvec[3]*hvec[3]) * invHmag3,
                                -(hvec[2]*hvec[3]) * invHmag3,
                                -(hvec[1]*hvec[3]) * invHmag3,
                                -(hvec[2]*hvec[3]) * invHmag3,
                                (hvec[1]*hvec[1] + hvec[2]*hvec[2]) * invHmag3)

        s3  = SMatrix{3,6}(-dhhatdh[1,2]*cart[6] + dhhatdh[1,3]*cart[5],
                           -dhhatdh[2,2]*cart[6] + dhhatdh[2,3]*cart[5],
                           -dhhatdh[3,2]*cart[6] + dhhatdh[3,3]*cart[5],
                            dhhatdh[1,1]*cart[6] - dhhatdh[1,3]*cart[4],
                            dhhatdh[2,1]*cart[6] - dhhatdh[2,3]*cart[4],
                            dhhatdh[3,1]*cart[6] - dhhatdh[3,3]*cart[4],
                           -dhhatdh[1,1]*cart[5] + dhhatdh[1,2]*cart[4],
                           -dhhatdh[2,1]*cart[5] + dhhatdh[2,2]*cart[4],
                           -dhhatdh[3,1]*cart[5] + dhhatdh[3,2]*cart[4],
                            dhhatdh[1,2]*cart[3] - dhhatdh[1,3]*cart[2],
                            dhhatdh[2,2]*cart[3] - dhhatdh[2,3]*cart[2],
                            dhhatdh[3,2]*cart[3] - dhhatdh[3,3]*cart[2],
                           -dhhatdh[1,1]*cart[3] + dhhatdh[1,3]*cart[1],
                           -dhhatdh[2,1]*cart[3] + dhhatdh[2,3]*cart[1],
                           -dhhatdh[3,1]*cart[3] + dhhatdh[3,3]*cart[1],
                            dhhatdh[1,1]*cart[2] - dhhatdh[1,2]*cart[1],
                            dhhatdh[2,1]*cart[2] - dhhatdh[2,2]*cart[1],
                            dhhatdh[3,1]*cart[2] - dhhatdh[3,2]*cart[1])
    end

    # Compute slice 2 (∂A[:,2] / ∂cart or ∂ithat / ∂cart)
    s2 = begin
        dtdr    = SMatrix{3,3}(-cart[2]*cart[5] - cart[3]*cart[6],
                                hvec[3] + cart[1]*cart[5],
                                cart[1]*cart[6] - hvec[2],
                                cart[2]*cart[4] - hvec[3],
                            -cart[1]*cart[4] - cart[3]*cart[6],
                                hvec[1] + cart[2]*cart[6],
                                hvec[2] + cart[3]*cart[4],
                                cart[3]*cart[5] - hvec[1],
                            -cart[1]*cart[4] - cart[2]*cart[5]) 

        dtdv    = SMatrix{3,3}(cart[2]*cart[2] + cart[3]*cart[3],
                            -cart[1]*cart[2],
                            -cart[1]*cart[3],
                            -cart[1]*cart[2],
                            cart[1]*cart[1] + cart[3]*cart[3],
                            -cart[2]*cart[3],
                            -cart[1]*cart[3],
                            -cart[2]*cart[3],
                            cart[1]*cart[1] + cart[2]*cart[2]) 

        tvec        = SVector(hvec[2]*cart[3] - hvec[3]*cart[2],
                            hvec[3]*cart[1] - hvec[1]*cart[3],
                            hvec[1]*cart[2] - hvec[2]*cart[1])
        tmag        = sqrt(tvec[1]*tvec[1] + tvec[2]*tvec[2] + tvec[3]*tvec[3])
        invTmag     = 1.0 / tmag
        invTmag3    = invTmag * invTmag * invTmag

        dthatdt = SMatrix{3,3}(invTmag - tvec[1]*tvec[1]*invTmag3,
                                -tvec[1]*tvec[2]*invTmag3,
                                -tvec[1]*tvec[3]*invTmag3,
                                -tvec[1]*tvec[2]*invTmag3,
                                invTmag - tvec[2]*tvec[2]*invTmag3,
                                -tvec[2]*tvec[3]*invTmag3,
                                -tvec[1]*tvec[3]*invTmag3,
                                -tvec[2]*tvec[3]*invTmag3,
                                invTmag - tvec[3]*tvec[3]*invTmag3)

        s2  = SMatrix{3,6}(dthatdt[1,1]*dtdr[1,1] + dthatdt[1,2]*dtdr[2,1] + dthatdt[1,3]*dtdr[3,1],
                        dthatdt[2,1]*dtdr[1,1] + dthatdt[2,2]*dtdr[2,1] + dthatdt[2,3]*dtdr[3,1],
                        dthatdt[3,1]*dtdr[1,1] + dthatdt[3,2]*dtdr[2,1] + dthatdt[3,3]*dtdr[3,1],
                        dthatdt[1,1]*dtdr[1,2] + dthatdt[1,2]*dtdr[2,2] + dthatdt[1,3]*dtdr[3,2],
                        dthatdt[2,1]*dtdr[1,2] + dthatdt[2,2]*dtdr[2,2] + dthatdt[2,3]*dtdr[3,2],
                        dthatdt[3,1]*dtdr[1,2] + dthatdt[3,2]*dtdr[2,2] + dthatdt[3,3]*dtdr[3,2],
                        dthatdt[1,1]*dtdr[1,3] + dthatdt[1,2]*dtdr[2,3] + dthatdt[1,3]*dtdr[3,3],
                        dthatdt[2,1]*dtdr[1,3] + dthatdt[2,2]*dtdr[2,3] + dthatdt[2,3]*dtdr[3,3],
                        dthatdt[3,1]*dtdr[1,3] + dthatdt[3,2]*dtdr[2,3] + dthatdt[3,3]*dtdr[3,3],
                        dthatdt[1,1]*dtdv[1,1] + dthatdt[1,2]*dtdv[2,1] + dthatdt[1,3]*dtdv[3,1],
                        dthatdt[2,1]*dtdv[1,1] + dthatdt[2,2]*dtdv[2,1] + dthatdt[2,3]*dtdv[3,1],
                        dthatdt[3,1]*dtdv[1,1] + dthatdt[3,2]*dtdv[2,1] + dthatdt[3,3]*dtdv[3,1],
                        dthatdt[1,1]*dtdv[1,2] + dthatdt[1,2]*dtdv[2,2] + dthatdt[1,3]*dtdv[3,2],
                        dthatdt[2,1]*dtdv[1,2] + dthatdt[2,2]*dtdv[2,2] + dthatdt[2,3]*dtdv[3,2],
                        dthatdt[3,1]*dtdv[1,2] + dthatdt[3,2]*dtdv[2,2] + dthatdt[3,3]*dtdv[3,2],
                        dthatdt[1,1]*dtdv[1,3] + dthatdt[1,2]*dtdv[2,3] + dthatdt[1,3]*dtdv[3,3],
                        dthatdt[2,1]*dtdv[1,3] + dthatdt[2,2]*dtdv[2,3] + dthatdt[2,3]*dtdv[3,3],
                        dthatdt[3,1]*dtdv[1,3] + dthatdt[3,2]*dtdv[2,3] + dthatdt[3,3]*dtdv[3,3])
    end
    
    return (s1,s2,s3)
end

function temp_computeA3bRTN(cart,t,p)
    # Compute rotation from MEE to Inertial
    mee2I   = mee2Inertial(cart)

    # Get thirdy body perterbations in inertial frame
    atbI    = computeThirdBodyPerterbations(cart, t, p)

    # Rotate perterbations to RTN frame
    Δr      = mee2I[1,1]*atbI[1] + mee2I[2,1]*atbI[2] + mee2I[3,1]*atbI[3]
    Δt      = mee2I[1,2]*atbI[1] + mee2I[2,2]*atbI[2] + mee2I[3,2]*atbI[3]
    Δn      = mee2I[1,3]*atbI[1] + mee2I[2,3]*atbI[2] + mee2I[3,3]*atbI[3]

    return SVector(Δr,Δt,Δn)
end

function computePartialOfCostIntegrandWrtMEEState(mee, p::MEEParams{MinimumFuel,T}) where {T}
    return SVector{7}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function computePartialOfCostIntegrandWrtMEEState(mee, p::MEEParams{MinimumFuelMayer,T}) where {T}
    return SVector{7}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function computePartialOfCostIntegrandWrtMEEState(mee, p::MEEParams{MinimumEnergyToFuel,T}) where {T}
    return SVector{7}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function computePartialOfCostIntegrandWrtMEEState(mee, p::MEEParams{MinimumTime,T}) where {T}
    return SVector{7}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end