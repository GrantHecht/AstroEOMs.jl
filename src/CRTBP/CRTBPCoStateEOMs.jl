# CRTBP CoState Dynamics
function crtbpMinimumFuelCoStateEOMs(u::AbstractArray, p::Tuple{SimpleSpacecraft,CRTBPParams}, γ)
    @inbounds begin
        # Get requirements 
        TU  = p[2].TU
        LU  = p[2].LU
        MU  = p[2].MU 
        μ   = p[2].μ

        # Scale requirements
        tMaxSc = p[1].tMax * TU * TU / (MU*LU*1000.0)

        # Compute Requirements
        xpmu    = u[1] + μ
        xpmum1  = u[1] + μ - 1
        r1      = sqrt(xpmu*xpmu + u[2]*u[2] + u[3]*u[3])
        r2      = sqrt(xpmum1*xpmum1 + u[2]*u[2] + u[3]*u[3])
        invr13  = r1^(-3)
        invr15  = r1^(-5)
        invr23  = r2^(-3)
        invr25  = r2^(-5)
        λv      = sqrt(u[11]^2 + u[12]^2 + u[13]^2)

        # Compute G 
        G11     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*(u[1] + μ)^2*invr15 - 
                    μ*invr23 + 3*μ*(u[1] + μ - 1)^2*invr25;
        G22     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*u[2]^2*invr15 - 
                    μ*invr23 + 3*μ*u[2]^2*invr25;
        G33     = -(1 - μ)*invr13 + 3*(1 - μ)*u[3]^2*invr15 - 
                    μ*invr23 + 3*μ*u[3]^2*invr25;
        G12     = 3*(1 - μ)*(u[1] + μ)*u[2]*invr15 + 
                    3*μ*(u[1] + μ - 1)*u[2]*invr25;
        G13     = 3*(1 - μ)*(u[1] + μ)*u[3]*invr15 + 
                    3*μ*(u[1] + μ - 1)*u[3]*invr25;
        G23     = 3*(1 - μ)*u[2]*u[3]*invr15 + 
                    3*μ*u[2]*u[3]*invr25;

        # Compute and return Dynamics
        dλ = @SVector [ -G11*u[11] - G12*u[12] - G13*u[13],
                        -G12*u[11] - G22*u[12] - G23*u[13],
                        -G13*u[11] - G23*u[12] - G33*u[13],
                        -u[8] + 2.0*u[12],
                        -u[9] - 2.0*u[11],
                        -u[10],
                        -λv*γ*tMaxSc / (u[7]*u[7])]
    end

    return dλ
end

function crtbpMinimumFuelCoStateEOMs!(dλ::AbstractArray, u::AbstractArray, p::Tuple{SimpleSpacecraft,CRTBPParams}, γ)
    @inbounds begin
    # Get requirements 
    TU  = p[2].TU
    LU  = p[2].LU
    MU  = p[2].MU 
    μ   = p[2].μ

    # Scale requirements
    tMaxSc = p[1].tMax * TU^2 / (MU*LU*1000.0)

    # Compute Requirements
    r1      = sqrt((u[1] + μ)^2 + u[2]^2 + u[3]^2)
    r2      = sqrt((u[1] + μ - 1)^2 + u[2]^2 + u[3]^2)
    invr13  = r1^(-3)
    invr15  = r1^(-5)
    invr23  = r2^(-3)
    invr25  = r2^(-5)
    λv      = sqrt(u[11]^2 + u[12]^2 + u[13]^2)

    # Compute G 
    G11     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*(u[1] + μ)^2*invr15 - 
                μ*invr23 + 3*μ*(u[1] + μ - 1)^2*invr25;
    G22     = 1 - (1 - μ)*invr13 + 3*(1 - μ)*u[2]^2*invr15 - 
                μ*invr23 + 3*μ*u[2]^2*invr25;
    G33     = -(1 - μ)*invr13 + 3*(1 - μ)*u[3]^2*invr15 - 
                μ*invr23 + 3*μ*u[3]^2*invr25;
    G12     = 3*(1 - μ)*(u[1] + μ)*u[2]*invr15 + 
                3*μ*(u[1] + μ - 1)*u[2]*invr25;
    G13     = 3*(1 - μ)*(u[1] + μ)*u[3]*invr15 + 
                3*μ*(u[1] + μ - 1)*u[3]*invr25;
    G23     = 3*(1 - μ)*u[2]*u[3]*invr15 + 
                3*μ*u[2]*u[3]*invr25;

    # Compute and return Dynamics
    dλ[1] = -G11*u[11] - G12*u[12] - G13*u[13]
    dλ[2] = -G12*u[11] - G22*u[12] - G23*u[13]
    dλ[3] = -G13*u[11] - G23*u[12] - G33*u[13]
    dλ[4] = -u[8] + 2.0*u[12]
    dλ[5] = -u[9] - 2.0*u[11]
    dλ[6] = -u[10]
    dλ[7] = -λv*γ*tMaxSc / (u[7]^2)
    end

    return @SVector [G11, G22, G33, G12, G13, G23]
end

