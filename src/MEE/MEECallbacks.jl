
function meePositivePCondition(x, t, integrator)
    return x[1] <= 0.0
end

function meePositivePAffect!(integrator)
    terminate!(integrator)
end

function meeMinimumFuelSwitchingCallbackCondition(x, t, integrator)

end