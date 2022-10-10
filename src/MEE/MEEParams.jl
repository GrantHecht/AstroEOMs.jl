# Modified Equinoctual Element EOM Parameters
struct MEEParams{costFlag <: CostFlag,tbdt} <: AbstractEOMParams
    # Initial epoch (TDB)
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
end

function MEEParams(initEpoch; μ = 3.986004415e5, primaryBodyID = 399, ref = "J2000",
                    thirdBodyEphemerides = nothing, 
                    costFunction::Type = MinimumFuel)

    # Check if we'll be using any perterbations
    perterbations = false
    if thirdBodyEphemerides !== nothing # Or other perts (use to set if any perts are on)
        perterbations = true
    end

    # Set flag to indicate if using 3rd body perterbations
    thirdBodyPerterbations = thirdBodyEphemerides !== nothing

    return MEEParams{costFunction, typeof(thirdBodyEphemerides)}(initEpoch, μ, primaryBodyID, ref, 
        perterbations, thirdBodyPerterbations, thirdBodyEphemerides)
end