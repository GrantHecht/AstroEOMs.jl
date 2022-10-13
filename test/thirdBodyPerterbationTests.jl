
using AstroEOMs, AstroUtils, ForwardDiff, SPICE, StaticArrays

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

# Define cartesian state
cart0       = SVector(-4.482091357854554e3,5.084589601512431e3,0.0,
                    -0.007552602581471e3,-0.006657701803705e3,0.0,1500.0)

# ===== Compute partials of qk
s               = getPosition(meeParams.thirdBodyEphemerides.ephems[1], meeParams.initEpoch)
dqkdr           = AstroEOMs.computePartialQkWrtR(view(cart0, 1:3), s)
dqkdrFD         = ForwardDiff.gradient(x -> AstroEOMs.computeQk(x,s), view(cart0, 1:3))
#display((dqkdr .- dqkdrFD) ./ eps.(dqkdr))
for i in eachindex(dqkdr)
    @test (dqkdr[i] - dqkdrFD[i]) / eps(dqkdr[i]) < 10.0
end

# ===== Compute partials of F(qk)
dFkdr           = AstroEOMs.computePartialFkWrtR(view(cart0, 1:3), s)
dFkdrFD         = ForwardDiff.gradient(x -> AstroEOMs.computeFk(x,s), view(cart0, 1:3))
#display((dFkdr .- dFkdrFD) ./ eps.(dFkdr))
for i in eachindex(dFkdr)
    @test (dFkdr[i] - dFkdrFD[i]) / eps(dFkdr[i]) < 10.0
end

# ===== Compute partials of inertial acceleration
ab, dabdcart    = AstroEOMs.computeThirdBodyPerterbationsAndPartials(cart0, 0.0, meeParams)
dabdcartFD      = ForwardDiff.jacobian(x -> AstroEOMs.computeThirdBodyPerterbations(x, 0.0, meeParams), cart0)
#display((dabdcart .- dabdcartFD) ./ eps.(dabdcart))
for i in eachindex(dabdcart)
    @test (dabdcart[i] - dabdcartFD[i]) / eps(dabdcart[i]) < 10.0
end