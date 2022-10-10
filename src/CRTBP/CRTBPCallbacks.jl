
function crtbpMEMFSwitchingCallbackCondition(x, t, integrator)
    # Get Requirements
    TU  = integrator.p[2].TU
    LU  = integrator.p[2].LU
    isp = integrator.p[1].isp

    ϵ       = integrator.p[3].ϵ
    utype   = integrator.p[3].state

    # Compute Scaling
    cSc = isp*9.81*TU / (LU * 1000.0)
    
    # Compute required
    λv = sqrt(x[11]^2 + x[12]^2 + x[13]^2)
    S  = computeS(x, λv, cSc)

    # Switching condition 
    out = 0.0
    if ϵ != 0.0
        if  utype == 0
            out = S - ϵ
        elseif utype == 2
            out = S + ϵ
        else
            out = abs(S) - ϵ
        end
    else
        out = S
    end

    return out
end

function crtbpMEMFSwitchingCallbackAffect!(integrator)
    # Get Requirements
    TU  = integrator.p[2].TU
    LU  = integrator.p[2].LU
    isp = integrator.p[1].isp

    ϵ   = integrator.p[3].ϵ
    utype = integrator.p[3].state

    # Switching affect 
    if ϵ != 0.0 
        if utype != 1
            integrator.p[3].state = 1
        else
            cSc = isp*9.81*TU / (LU*1000.0)
            λv = sqrt(integrator.u[11]*integrator.u[11] + 
                      integrator.u[12]*integrator.u[12] + 
                      integrator.u[13]*integrator.u[13])
            S = computeS(integrator.u, λv, cSc)
            if S < 0.0
                integrator.p[3].state = 2
            else
                integrator.p[3].state = 0
            end
        end

    else 
        if utype == 0
            integrator.p[3].state = 2
        else
            integrator.p[3].state = 0
        end
    end
    return nothing
end