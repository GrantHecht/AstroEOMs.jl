using AstroEOMs, AstroUtils, ForwardDiff, SPICE, StaticArrays, LinearAlgebra

# Furnsh default SPICE kernels
furnshDefaults()

# Compute initial epoch
initEpoch   = utc2et("2022-10-07T12:00:00") 

# Construct ephemerides for third body perterbations
ephemTspan  = (initEpoch - 3600, initEpoch + 6.0*24.0*3600.0)
nPoints     = ceil(Int64, (ephemTspan[2] - ephemTspan[1]) / 1800.0)
targIDs     = [10,301]
tbEphems    = Ephemerides(ephemTspan, nPoints, targIDs, 399, "J2000")

# Construct Parameters
meeParams   = MEEParams(initEpoch; thirdBodyEphemerides = tbEphems)

# Create simple spacecraft
sp          = SimpleSpacecraft(1500.0, 1500.0, 10.0, 300.0)

# Define cartesian state
cart = SVector(6524.834, 6862.875,  6448.296, # [km] 
             4.901327, 5.533756, -1.976341) # [km/s]
mee0       = convertState(cart, Cartesian, MEE, meeParams.mu)
meeWithM   = SVector(mee0[1], mee0[2], mee0[3], mee0[4],
                mee0[5], mee0[6], 1500.0)

# ===== Test MEE perterbation partials
ab, dapcart    = AstroEOMs.meeComputePerterbationPartials(meeWithM, meeParams, 0.0)
dapcartFD      = ForwardDiff.jacobian(x -> AstroEOMs.meeComputePerterbations(x, meeParams, 0.0), meeWithM)
diff           = dapcart .- dapcartFD
for i in eachindex(diff)
    @test diff[i] / eps(dapcart[i]) < 75.0
end

# ===== Test MEE dynamics partials
α       = [1.0, 2.0, 3.0]
α      ./ norm(α)
tvec    = 10.0e-3 * 1.0 * α
dfdx    = AstroEOMs.meeComputeDynamicsPartials(meeWithM, meeParams, 0.0, tvec ./ meeWithM[7])
dfdxFD  = ForwardDiff.jacobian(x -> AstroEOMs.meeEomControl(x, meeParams, 0.0, tvec ./ x[7]), meeWithM)
diff    = (dfdx[1:6,:] .- dfdxFD) ./ eps.(dfdx[1:6,:])
for i in eachindex(diff)
    @test diff[i]  < 50.0
end

# ===== Test MEE co-state dynamics
function meeHamiltonian(u,ps,t,tvec)
    au  = tvec ./ u[7]
    fx  = AstroEOMs.meeEomControl(view(u,1:6), ps, t, au)
    fm  = -norm(tvec) / (300.0 * 9.81)
    return dot(view(u,8:13),view(fx,1:6)) + fm*u[14] - fm
end
y0      = SVector{14}(meeWithM[1], meeWithM[2], meeWithM[3], meeWithM[4], meeWithM[5], meeWithM[6], meeWithM[7],
                  10.0, 0.5, 0.25, 0.1, 0.22, 0.1, 1.1)
fλ      = AstroEOMs.meeCoStateEOMs(y0, meeParams, 0.0, tvec ./ y0[7])
fλFD    = -ForwardDiff.gradient(x -> meeHamiltonian(
                SVector{14}(x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                            10.0, 0.5, 0.25, 0.1, 0.22, 0.1, 1.1),
                meeParams, 0.0, tvec), meeWithM)
diff    = (fλ .- fλFD) ./ eps.(fλ)
for i in eachindex(diff)
    @test diff[i] < 20.0
end