
function meeMEMFInitializeSwitchingState!(sw, u, t, sp, mp)
    # Compute optimal thrust direction
    α, nBtλ = meeComputeOptimalThrustUnitVector(u, mp)

    # Compute switching
    S = meeComputeS(u, nBtλ, sp, mp)

    # Set switching state
    if S < -sw.ϵ
        sw.state = 0
    elseif S > sw.ϵ
        sw.state = 2
    else 
        sw.state = 1
    end
    return nothing
end

function meeMEMFSwitchingCallbackCondition(u, t, integrator)
    # Grab parameters
    sp = integrator.p[1]
    mp = integrator.p[2]
    sw = integrator.p[3]

    # Compute optimal thrust direction
    α, nBtλ = meeComputeOptimalThrustUnitVector(u, mp)

    # Compute switching
    S = meeComputeS(u, nBtλ, sp, mp)
    
    # Switching condition
    out = 0.0
    if sw.ϵ != 0.0
        if sw.state == 0
            out = S + sw.ϵ
        elseif sw.state == 2
            out = S - sw.ϵ
        else
            out = abs(S) - sw.ϵ
        end
    else
        out = S
    end

    return out
end

function meeMEMFSwitchingCallbackAffect!(integrator)
    # Grab parameters
    sp = integrator.p[1]
    mp = integrator.p[2]
    sw = integrator.p[3]

    # Switching affect
    if sw.ϵ != 0.0
        if sw.state != 1
            sw.state = 1
        else
            # Compute optimal thrust direction
            α, nBtλ = meeComputeOptimalThrustUnitVector(integrator.u, mp)

            # Compute switching
            S = meeComputeS(integrator.u, nBtλ, sp, mp)
            if S < 0.0
                sw.state = 0
            else
                sw.state = 2
            end
        end
    else
        if sw.state == 0
            sw.state = 2
        else
            sw.state = 0
        end
    end
    return nothing
end