module AstroEOMs

using StaticArrays
using ForwardDiff
using LinearAlgebra

using Pkg
Pkg.add(url="https://github.com/GrantHecht/AstroUtils.jl.git")
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
