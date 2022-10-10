# A simple struct to be used for handling switching with 
# optimal control variational equation integration
mutable struct SwitchingStruct 
    # Continuation parameter
    ϵ::Float64

    # Switch state flag
    state::Int
end

SwitchingStruct(state) = SwitchingStruct(0.0, state)