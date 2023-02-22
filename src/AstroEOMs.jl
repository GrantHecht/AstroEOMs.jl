module AstroEOMs

using StaticArrays
using ForwardDiff
using LinearAlgebra
using AstroUtils

include("./AbstractEOMParams.jl")
include("./SimpleSpacecraft.jl")
include("./CostFlags.jl")
include("./SwitchingStruct.jl")

# Include force models
include("./ForceModels/ForceModels.jl")

# CRTBP Imports
include("./CRTBP/CRTBP.jl")

# MEE Imports 
include("./MEE/MEE.jl")

export SimpleSpacecraft, SwitchingStruct

end
