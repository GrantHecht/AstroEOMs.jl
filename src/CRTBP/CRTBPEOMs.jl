# CR3BP Eq. of Motion w/o control out-of-place with SVector
function crtbpEomNoControl(u::AbstractArray, p::CRTBPParams, t)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    return @SVector [u[4], u[5], u[6],
                     gx + hx,
                     gy + hy,
                     gz]
end

# CRTBP Eq. of Motion w/o control in-place
function crtbpEomNoControl!(du::AbstractArray,u::AbstractArray,p::CRTBPParams,t)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Set du 
    du[1:3] .= view(u,4:6)
    du[4] = gx + hx 
    du[5] = gy + hy
    du[6] = gz
    return nothing
end

# CR3BP Eq. of Motion w/ control out-of-place with SVector
function crtbpEomControl(u::AbstractArray,p::CRTBPParams,t,a::AbstractArray)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = inv(r1^3)
    invR23 = inv(r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Compute and return dynamics
    dx  = u[4]
    dy  = u[5]
    dz  = u[6]
    ddx = gx + hx + a[1]
    ddy = gy + hy + a[2]
    ddz = gz + a[3]

    return @SVector [dx,dy,dz,ddx,ddy,ddz]
end

# CRTBP Eq. of Motion w/ control in-place 
function crtbpEomControl!(du::AbstractArray,u::AbstractArray,p::CRTBPParams,t,a::AbstractArray)
    # Precompute required
    xpμ = u[1] + p.μ
    r1  = sqrt(xpμ^2 + u[2]^2 + u[3]^2)
    r2  = sqrt((xpμ - 1)^2 + u[2]^2 + u[3]^2)

    # Compute g(r)
    invR13 = 1.0 / (r1^3)
    invR23 = 1.0 / (r2^3)
    gx = u[1] - (1.0 - p.μ)*xpμ*invR13 - p.μ*(xpμ - 1.0)*invR23
    gy = u[2] - (1.0 - p.μ)*u[2]*invR13 - p.μ*u[2]*invR23 
    gz = -(1.0 - p.μ)*u[3]*invR13 - p.μ*u[3]*invR23

    # Compute h(v)
    hx = 2.0*u[5]
    hy = -2.0*u[4]

    # Set du 
    du[1:3] .= view(u, 4:6)
    du[4] = gx + hx + a[1]
    du[5] = gy + hy + a[2]
    du[6] = gz + a[3]
end