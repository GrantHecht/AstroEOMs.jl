
function meeOptimalControlEOMs(u::AbstractVector, p::Tuple{SimpleSpacecraft,MEEParams,SwitchingStruct}, t)
    # Grab parameter structs
    sp = p[1]
    mp = p[2]
    sw = p[3]

    # Compute optimal thrust direction
    α, nBtλ = meeComputeOptimalThrustUnitVector(u, mp)

    # Compute throttling
    S  = meeComputeS(u, nBtλ, sp, mp)
    ut = meeComputeU(S, mp, sw)

    # Compute acceleration due to thrust
    auMag_us    = ut*sp.tMax / (1000.0 * mp.MU * u[7])
    au_us       = SVector(auMag_us*α[1], auMag_us*α[2], auMag_us*α[3])
    au          = scaleAcceleration(mp, au_us)

    # Compute derivaties
    dx = meeEomControl(view(u, 1:6), mp, t, au)
    dm = -ut*sp.tMax*mp.TU / (sp.c * mp.MU)
    dλ = meeCoStateEOMs(u, mp, t, au)

    return SVector{14}(dx[1], dx[2], dx[3], dx[4], dx[5], dx[6], dm,
                dλ[1], dλ[2], dλ[3], dλ[4], dλ[5], dλ[6], dλ[7])
end

function meeOptimalControlHamiltonian(u::AbstractVector, p::Tuple{SimpleSpacecraft,MEEParams,SwitchingStruct}, t)
    # Grab parameter structs
    sp = p[1]
    mp = p[2]
    sw = p[3]

    # Compute optimal thrust direction
    α, nBtλ = meeComputeOptimalThrustUnitVector(u, mp)

    # Compute throttling
    S  = meeComputeS(u, nBtλ, sp, mp)
    ut = meeComputeU(S, mp, sw)

    # Compute acceleration due to thrust
    auMag = ut*sp.tMax*mp.TU*mp.TU / (1000.0 * mp.LU * mp.MU * u[7])
    au = SVector(auMag*α[1], auMag*α[2], auMag*α[3])

    # Compute derivaties
    dx = meeEomControl(view(u, 1:6), mp, t, au)
    dm = -ut*sp.tMax*mp.TU / (sp.c * mp.MU)
    #dλ = meeCoStateEOMs(u, mp, t, au)

    return dx[1]*u[8] + dx[2]*u[9] + dx[3]*u[10] + dx[4]*u[11] + dx[5]*u[12] + dx[6]*u[13] + dm*u[14] + meeComputeCostIntegrand(ut,mp,sw)
end

function meeComputeOptimalThrustUnitVector(u::AbstractVector, p::MEEParams)
    # Get scaled gravitational parameter
    mu = getScaledGravityParameter(p)

    # Compute requirements
    w  = 1.0 + u[2]*cos(u[6]) + u[3]*sin(u[6])
    ss = 1.0 + u[4]*u[4] + u[5]*u[5]
    if u[1] <= 0.0
        display(u[1])
    end
    sqrtPmu = sqrt(u[1] / mu)
    wInv    = 1.0 / w

    # Compute B matrix
    B  = SMatrix{6, 3}(0.0, sqrtPmu*sin(u[6]), -sqrtPmu*cos(u[6]), 0.0, 0.0, 0.0,
                       2.0*u[1]*wInv*sqrtPmu, sqrtPmu*wInv*((w + 1.0)*cos(u[6]) + u[2]), sqrtPmu*wInv*((w + 1.0)*sin(u[6]) + u[3]), 0.0, 0.0, 0.0,
                       0.0, -sqrtPmu*u[3]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])), sqrtPmu*u[2]*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])),
                       0.5*sqrtPmu*ss*wInv*cos(u[6]), 0.5*sqrtPmu*ss*wInv*sin(u[6]), sqrtPmu*wInv*(u[4]*sin(u[6]) - u[5]*cos(u[6])))

    # Compute Bᵀλₓ
    Btλ = SVector(B[1,1]*u[8] + B[2,1]*u[9] + B[3,1]*u[10] + B[4,1]*u[11] + B[5,1]*u[12] * B[6,1]*u[13],
                  B[1,2]*u[8] + B[2,2]*u[9] + B[3,2]*u[10] + B[4,2]*u[11] + B[5,2]*u[12] * B[6,2]*u[13],
                  B[1,3]*u[8] + B[2,3]*u[9] + B[3,3]*u[10] + B[4,3]*u[11] + B[5,3]*u[12] * B[6,3]*u[13])

    # Compute norm of Btλ
    nBtλ = sqrt(Btλ[1]*Btλ[1] + Btλ[2]*Btλ[2] + Btλ[3]*Btλ[3])

    # Return optimal thrust direction and nBtλ
    return (SVector(-Btλ[1] / nBtλ, -Btλ[2] / nBtλ, -Btλ[3] / nBtλ), nBtλ)
end

function meeComputeCostIntegrand(ut, mp::MEEParams{MinimumTime,T},sw::SwitchingStruct) where {T}
    return 1.0
end

function meeComputeCostIntegrand(ut, mp::MEEParams{MinimumFuel,T},sw::SwitchingStruct) where {T}
    return ut
end

function meeComputeCostIntegrand(ut, mp::MEEParams{MinimumEnergyToFuel,T},sw::SwitchingStruct) where {T}
    return ut - sw.ϵ*ut*(1.0 - ut) 
end

function meeComputeCostIntegrand(ut, mp::MEEParams{MinimumFuelMayer,T},sw::SwitchingStruct) where {T}
    return 0.0
end

function meeComputeS(u, nBtλ, sp::SimpleSpacecraft, mp::MEEParams{MinimumTime, T}) where {T}
    c_sc = sp.c * mp.TU / (mp.LU * 1000.0)
    return c_sc*nBtλ / u[7] + u[14] 
end

function meeComputeS(u, nBtλ, sp::SimpleSpacecraft, mp::MEEParams{MinimumFuel, T}) where {T}
    c_sc = sp.c * mp.TU / (mp.LU * 1000.0)
    return c_sc * nBtλ / u[7] + u[14]  - 1.0
end

function meeComputeS(u, nBtλ, sp::SimpleSpacecraft, mp::MEEParams{MinimumFuelMayer, T}) where {T}
    c_sc = sp.c * mp.TU / (mp.LU * 1000.0)
    return c_sc * nBtλ / u[7] + u[14] 
end

function meeComputeS(u, nBtλ, sp::SimpleSpacecraft, mp::MEEParams{MinimumEnergyToFuel, T}) where {T}
    c_sc = sp.c * mp.TU / (mp.LU * 1000.0)
    return c_sc * nBtλ / u[7] + u[14]  - 1.0
end

function meeComputeU(S, mp::MEEParams{MinimumTime, T}, sw::SwitchingStruct) where {T}
    return 0.5*(1.0 + tanh(S / sw.ϵ))
end

function meeComputeU(S, mp::MEEParams{MinimumFuel, T}, sw::SwitchingStruct) where {T}
    return 0.5*(1.0 + tanh(S / sw.ϵ))
end

function meeComputeU(S, mp::MEEParams{MinimumFuelMayer, T}, sw::SwitchingStruct) where {T}
    return 0.5*(1.0 + tanh(S / sw.ϵ))
end

function meeComputeU(S, mp::MEEParams{MinimumEnergyToFuel, T}, sw::SwitchingStruct) where {T}
    if sw.state == 0
        return 0.0
    elseif sw.state == 2
        return 1.0
    else
        return (S + sw.ϵ) / (2.0*sw.ϵ)
    end
end