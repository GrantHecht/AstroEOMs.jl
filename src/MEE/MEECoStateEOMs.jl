
function meeCoStateEOMs(u::AbstractArray, p::MEEParams, t, au) 
    # Compute partials of MEE EOMs w.r.t. MEE state
    dfdx = meeComputeDynamicsStatePartials(u, p, t, au)

    # Compute partials of Integral Cost Term w.r.t. MEE state
    dLdx = computePartialOfCostIntegrandWrtMEEState(u, p)

    # Compute and return co-state dynamics
    return SVector{7}(-u[8]*dfdx[1,1] - u[9]*dfdx[2,1] - u[10]*dfdx[3,1] - u[11]*dfdx[4,1] - u[12]*dfdx[5,1] - u[13]*dfdx[6,1] - u[14]*dfdx[7,1] - dLdx[1],
                      -u[8]*dfdx[1,2] - u[9]*dfdx[2,2] - u[10]*dfdx[3,2] - u[11]*dfdx[4,2] - u[12]*dfdx[5,2] - u[13]*dfdx[6,2] - u[14]*dfdx[7,2] - dLdx[2],
                      -u[8]*dfdx[1,3] - u[9]*dfdx[2,3] - u[10]*dfdx[3,3] - u[11]*dfdx[4,3] - u[12]*dfdx[5,3] - u[13]*dfdx[6,3] - u[14]*dfdx[7,3] - dLdx[3],
                      -u[8]*dfdx[1,4] - u[9]*dfdx[2,4] - u[10]*dfdx[3,4] - u[11]*dfdx[4,4] - u[12]*dfdx[5,4] - u[13]*dfdx[6,4] - u[14]*dfdx[7,4] - dLdx[4],
                      -u[8]*dfdx[1,5] - u[9]*dfdx[2,5] - u[10]*dfdx[3,5] - u[11]*dfdx[4,5] - u[12]*dfdx[5,5] - u[13]*dfdx[6,5] - u[14]*dfdx[7,5] - dLdx[5],
                      -u[8]*dfdx[1,6] - u[9]*dfdx[2,6] - u[10]*dfdx[3,6] - u[11]*dfdx[4,6] - u[12]*dfdx[5,6] - u[13]*dfdx[6,6] - u[14]*dfdx[7,6] - dLdx[6],
                      -u[8]*dfdx[1,7] - u[9]*dfdx[2,7] - u[10]*dfdx[3,7] - u[11]*dfdx[4,7] - u[12]*dfdx[5,7] - u[13]*dfdx[6,7] - u[14]*dfdx[7,7] - dLdx[7])
end