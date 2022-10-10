
function meeCoStateEOMs(u::AbstractArray, p::Tuple{SimpleSpacecraft,MEEParams}, au) 
    # Compute perterbations and partials
    
end

function meeComputePerterbationPartials(u::AbstractArray, p::Tuple{SimpleSpacecraft,MEEParams})
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
        cart        = convertState(u, MEE, Cartesian, p.mu)
        dcartdmee   = stateConvertPartials()
    end
end