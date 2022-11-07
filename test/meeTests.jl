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
meeParams   = MEEParams(initEpoch; #thirdBodyEphemerides = tbEphems, 
                    LU = 384400.0, MU = 1.0, TU = 24.0*3600.0)

# Create simple spacecraft
sp          = SimpleSpacecraft(1500.0, 1500.0, 10.0, 300.0)

# Define cartesian state
cart = SVector(6524.834, 6862.875,  6448.296, # [km] 
               4.901327, 5.533756, -1.976341) # [km/s]
mee0       = convertState(cart, Cartesian, MEE, meeParams.mu)
meeWithM   = SVector(mee0[1], mee0[2], mee0[3], mee0[4],
                mee0[5], mee0[6], 1500.0)
meeWithM_s = AstroEOMs.scaleStateWithMass(meeParams, meeWithM)

# ===== Test that all forms of MEE equations return the same states
dxnc    = AstroEOMs.meeEomNoControl(mee0, meeParams, 0.0)
dxncIP  = zeros(6)
AstroEOMs.meeEomNoControl!(dxncIP, mee0, meeParams, 0.0)
dxc     = AstroEOMs.meeEomControl(mee0, meeParams, 0.0, zeros(3))
dxcIP   = zeros(6)
AstroEOMs.meeEomControl!(dxcIP, mee0, meeParams, 0.0, zeros(3))

for i in eachindex(dxnc)
    @test dxnc[i] == dxncIP[i]
    @test dxnc[i] == dxc[i]
    @test dxnc[i] == dxcIP[i]
end

# ===== Test MEE perterbation partials
ab, dapcart    = AstroEOMs.meeComputePerterbationStatePartials(meeWithM, meeParams, 0.0)
dapdt          = AstroEOMs.meeComputePerterbationTimePartials(meeWithM, meeParams, 0.1)
dapcartFD      = ForwardDiff.jacobian(x -> AstroEOMs.meeComputePerterbations(x, meeParams, 0.0), meeWithM)
dapdtFD        = ForwardDiff.jacobian(t -> AstroEOMs.meeComputePerterbations(meeWithM, meeParams, t[1]), [0.1])

diff           = abs.(dapcart .- dapcartFD)
for i in eachindex(diff)
    @test diff[i] / eps(dapcart[i]) < 75.0
end
diff           = zeros(3)
for i in 1:3; diff[i] = abs(dapdt[i] - dapdtFD[i,1]); end
for i in eachindex(diff)
    @test diff[i] / eps(dapdt[i]) < 50.0
end

# ===== Test MEE dynamics partials
α       = [1.0, 2.0, 3.0]
α      ./ norm(α)
tvec_s  = (10.0e-3 * 1.0 * meeParams.TU*meeParams.TU / (meeParams.MU*meeParams.LU)) * α
dfdx    = AstroEOMs.meeComputeDynamicsStatePartials(meeWithM_s, meeParams, 0.0, tvec_s ./ meeWithM_s[7])
dfdt    = AstroEOMs.meeComputeDynamicsTimePartials(meeWithM_s, meeParams, 0.1)
dfdsf   = AstroEOMs.meeComputeDynamicsThrustScaleFactorPartials(meeWithM_s, (sp,meeParams), 0.0, α)
dfdα    = AstroEOMs.meeComputeDynamicsThrustUnitVectorPartials(meeWithM_s, (sp, meeParams), 0.0, 1.0)
dfdxFD  = ForwardDiff.jacobian(x -> AstroEOMs.meeEomControl(x, meeParams, 0.0, tvec_s ./ x[7]), meeWithM_s)
dfdtFD  = ForwardDiff.jacobian(t -> AstroEOMs.meeEomControl(meeWithM_s, meeParams, t[1], tvec_s ./ meeWithM_s[7]), [0.1])
dfdsfFD = ForwardDiff.jacobian(sf -> AstroEOMs.meeEomControl(meeWithM_s, meeParams, 0.0, 
            (sf[1] * sp.tMax * meeParams.TU^2 / (1000.0 * meeParams.MU * meeParams.LU * meeWithM_s[7])).*α), [1.0])
dfdαFD  = ForwardDiff.jacobian(uv -> AstroEOMs.meeEomControl(meeWithM_s, meeParams, 0.0, 
            (1.0 * sp.tMax * meeParams.TU^2 / (1000.0 * meeParams.MU * meeParams.LU * meeWithM_s[7])).*uv), α)

# State partial tests
diff    = abs.(dfdx[1:6,:] .- dfdxFD) ./ eps.(dfdx[1:6,:])
for i in eachindex(diff)
    @test diff[i]  < 50.0
end

# Time partial tests
diff    = zeros(6)
for i in 1:6; diff[i] = abs(dfdt[i] - dfdtFD[i,1]) / eps(dfdt[i]); end
for i in eachindex(diff)
    @test diff[i] < 50.0 
end

# Thrust scale factor tests
diff    = zeros(6)
for i in 1:6; diff[i] = abs(dfdsf[i] - dfdsfFD[i,1]) / eps(dfdsf[i]); end
for i in eachindex(diff)
    @test diff[i] < 50.0
end

# Thrust unit vector tests
diff    = abs.(dfdα[1:6,:] .- dfdαFD) ./ eps.(dfdα[1:6,:]) 
for i in eachindex(diff)
    @test diff[i] < 50.0
end

# ===== Test MEE co-state dynamics
function meeHamiltonian(u,ps,t,tvec)
    au  = tvec ./ u[7]
    fx  = AstroEOMs.meeEomControl(view(u,1:6), ps, t, au)
    fm  = -norm(tvec) / (300.0 * 9.81e-3 * ps.TU / ps.LU)
    return dot(view(u,8:13),view(fx,1:6)) + fm*u[14] - fm
end
y0_s    = SVector{14}(meeWithM_s[1], meeWithM_s[2], meeWithM_s[3], meeWithM_s[4], meeWithM_s[5], meeWithM_s[6], meeWithM_s[7],
                  10.0, 0.5, 0.25, 0.1, 0.22, 0.1, 1.1)

fλ      = AstroEOMs.meeCoStateEOMs(y0_s, meeParams, 0.1, tvec_s ./ y0_s[7])
fλFD    = -ForwardDiff.gradient(
                x -> meeHamiltonian(SVector{14}(x[1], x[2], x[3], x[4], x[5], x[6], x[7], 10.0, 0.5, 0.25, 0.1, 0.22, 0.1, 1.1),meeParams, 0.1, tvec_s), 
                meeWithM_s)
diff    = abs.(fλ .- fλFD) ./ eps.(fλ)

for i in eachindex(diff)
    @test diff[i] < 10.0
end