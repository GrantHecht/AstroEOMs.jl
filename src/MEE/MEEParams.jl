# Modified Equinoctual Element EOM Parameters
struct MEEParams{costFlag <: CostFlag,tbdt} <: AbstractEOMParams
    # Initial epoch (TDB sec)
    initEpoch::Float64

    # Primary body gravitational parameter
    mu::Float64

    # Primary body NAIF ID
    primaryBodyID::Int

    # Reference frame
    ref::String

    # ==== Perterbation flags, settings, and data
    perterbations::Bool

    # == Third body perterbations
    thirdBodyPerterbations::Bool  

    # Third body perterbation data
    thirdBodyEphemerides::tbdt

    # == Nonspherical gravity
    nonsphericalGravity::Bool
    onlyJ2::Bool
    J2::Float64     # J2 Coefficnent
    Rcb::Float64    # Radius of central body

    # Problem units
    TU::Float64 # [sec]
    LU::Float64 # [km]
    MU::Float64 # [kg]
end

function MEEParams(initEpoch; μ = 3.986004415e5, primaryBodyID = 399, ref = "J2000",
                   TU = 1.0, LU = 1.0, MU = 1.0,
                   thirdBodyEphemerides = nothing, 
                   nonsphericalGravity = false,
                   J2 = 1086.639e-6,
                   Rcb = 6378.14,
                   costFunction::Type = MinimumFuel)

    # Check if we'll be using any perterbations
    perterbations = false
    if thirdBodyEphemerides !== nothing # Or other perts (use to set if any perts are on)
        perterbations = true
    end

    # Set flag to indicate if using 3rd body perterbations
    thirdBodyPerterbations = thirdBodyEphemerides !== nothing

    return MEEParams{costFunction, typeof(thirdBodyEphemerides)}(initEpoch, μ, primaryBodyID, ref, 
        perterbations, thirdBodyPerterbations, thirdBodyEphemerides, nonsphericalGravity, true, 
        J2, Rcb, TU, LU, MU)
end

# Get scaled gravitational parameter
function getScaledGravityParameter(mp::MEEParams)
    sf = mp.TU*mp.TU / (mp.LU*mp.LU*mp.LU)
    return mp.mu*sf 
end

# Scale time
scaleTime(mp::MEEParams, t) = t / mp.TU

# Unscale time
unscaleTime(mp::MEEParams, t) = t * mp.TU

# Scales acceleration with units km / s^2 
function scaleAcceleration(mp::MEEParams, a_unscaled)
    sf = mp.TU*mp.TU / mp.LU
    return SVector{3}(a_unscaled[1]*sf, a_unscaled[2]*sf, a_unscaled[3]*sf)
end

# Unscale acceleration
function unscaleAcceleration(mp::MEEParams, a_scaled)
    sf = mp.LU / (mp.TU * mp.TU)
    return SVector{3}(a_scaled[1]*sf, a_scaled[2]*sf, a_scaled[3]*sf)
end

# Unscale MEE state
function unscaleState(mp::MEEParams, state)
    return SVector{6}(mp.LU*state[1], state[2], state[3], state[4], state[5], state[6])
end

# Scale MEE state
function scaleState(mp::MEEParams, state)
    return SVector{6}(state[1] / mp.LU, state[2], state[3], state[4], state[5], state[6])
end

# Unscale MEE state with mass
function unscaleStateWithMass(mp::MEEParams, state)
    return SVector{7}(mp.LU*state[1], state[2], state[3], state[4], state[5], state[6], mp.MU*state[7])
end

# Scale MEE state with mass
function scaleStateWithMass(mp::MEEParams, state)
    return SVector{7}(state[1] / mp.LU, state[2], state[3], state[4], state[5], state[6], state[7] / mp.MU)
end

# Scale MEE dynamics partials w.r.t MEE state
function scaleMEEDynamicsPartials(mp::MEEParams, dfdx)
    dpsf    = mp.TU / mp.LU
    dmsf    = mp.TU / mp.MU
    dndpsf  = mp.TU * mp.LU
    dndmsf  = mp.TU * mp.MU
    dmdpsf  = mp.TU * mp.LU / mp.MU 
    dpdmsf  = mp.TU * mp.MU / mp.LU

    return SMatrix{7,7}(dfdx[1,1]*mp.TU,   dfdx[2,1]*dndpsf, dfdx[3,1]*dndpsf, dfdx[4,1]*dndpsf, dfdx[5,1]*dndpsf, dfdx[6,1]*dndpsf, dfdx[7,1]*dmdpsf,
                        dfdx[1,2]*dpsf,    dfdx[2,2]*mp.TU,  dfdx[3,2]*mp.TU,  dfdx[4,2]*mp.TU,  dfdx[5,2]*mp.TU,  dfdx[6,2]*mp.TU,  dfdx[7,2]*dmsf,
                        dfdx[1,3]*dpsf,    dfdx[2,3]*mp.TU,  dfdx[3,3]*mp.TU,  dfdx[4,3]*mp.TU,  dfdx[5,3]*mp.TU,  dfdx[6,3]*mp.TU,  dfdx[7,2]*dmsf,
                        dfdx[1,4]*dpsf,    dfdx[2,4]*mp.TU,  dfdx[3,4]*mp.TU,  dfdx[4,4]*mp.TU,  dfdx[5,4]*mp.TU,  dfdx[6,4]*mp.TU,  dfdx[7,2]*dmsf,
                        dfdx[1,5]*dpsf,    dfdx[2,5]*mp.TU,  dfdx[3,5]*mp.TU,  dfdx[4,5]*mp.TU,  dfdx[5,5]*mp.TU,  dfdx[6,5]*mp.TU,  dfdx[7,2]*dmsf,
                        dfdx[1,6]*dpsf,    dfdx[2,6]*mp.TU,  dfdx[3,6]*mp.TU,  dfdx[4,6]*mp.TU,  dfdx[5,6]*mp.TU,  dfdx[6,6]*mp.TU,  dfdx[7,2]*dmsf,
                        dfdx[1,7]*dpdmsf,  dfdx[2,7]*dndmsf, dfdx[3,7]*dndmsf, dfdx[4,7]*dndmsf, dfdx[5,7]*dndmsf, dfdx[6,7]*dndmsf, dfdx[7,2]*mp.TU)
end