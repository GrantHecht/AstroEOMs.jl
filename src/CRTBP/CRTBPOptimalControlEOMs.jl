# CR3BP Indirect EOMs out-of-place with SVector
function crtbpMinimumFuelOptimalControlMEMFEOMs(u::AbstractVector, p::Tuple{SimpleSpacecraft,CRTBPParams,SwitchingStruct}, t)
    @inbounds begin

    # Get requirements 
    TU  = p[2].TU
    LU  = p[2].LU
    VU  = p[2].VU 
    MU  = p[2].MU 
    isp = p[1].isp 
    c   = p[1].c

    # Scale requirements
    tMaxSc = p[1].tMax * TU * TU / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]*u[11] + u[12]*u[12] + u[13]*u[13])
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p[3].state, p[3].ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Derivatives
    dx = crtbpEomControl(u, p[2], t, at)
    dλ = crtbpMinimumFuelCoStateEOMs(u, (p[1], p[2]), γ)
    du = @SVector [dx[1], dx[2], dx[3], dx[4], dx[5], dx[6], -γ*tMaxSc / cSc,
                   dλ[1], dλ[2], dλ[3], dλ[4], dλ[5], dλ[6], dλ[7]]

    return du
    end
end

function crtbpMinimumFuelOptimalControlMEMFEOMs!(du::AbstractVector, u::AbstractVector, 
                                                p::Tuple{SimpleSpacecraft,CRTBPParams,SwitchingStruct}, t)
    @inbounds begin
    # Get requirements 
    TU  = p[2].TU
    LU  = p[2].LU
    VU  = p[2].VU 
    MU  = p[2].MU 
    isp = p[1].isp 
    c   = p[1].c

    # Scale requirements
    tMaxSc = p[1].tMax * TU^2 / (MU*LU*1000.0)
    cSc = c*TU / (LU*1000.0)

    # Compute thrust direction
    λv = sqrt(u[11]^2 + u[12]^2 + u[13]^2)
    invλv = 1.0 / λv 
    α = @SVector [-u[11]*invλv, -u[12]*invλv, -u[13]*invλv]

    # Compute throttline 
    S = computeS(u, λv, cSc)
    γ = computeU(S, p[3].state, p[3].ϵ)

    # Compute Thrust Acceleration
    atMag = γ*tMaxSc / u[7]
    at = @SVector [α[1]*atMag,
                   α[2]*atMag,
                   α[3]*atMag]

    # Compute dynamics
    crtbpEomControl!(view(du,1:6), u, p[2], t, at)
    du[7] = -γ*tMaxSc / cSc
    GVec  = crtbpMinimumFuelCoStateEOMs!(view(du,8:14), u, (p[1], p[2]), γ)
    end
    return GVec
end

# Utility Functions 
function computeS(x::AbstractVector, λv, cSc)
    @inbounds begin
        return -λv*cSc / x[7] - x[14] + 1.0
    end
end

function computeU(S, utype, ϵ)
    if utype == 2
        return 1.0
    elseif utype == 0
        return 0.0
    else
        return (ϵ - S) / (2.0*ϵ)
    end
end