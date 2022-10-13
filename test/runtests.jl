using AstroEOMs
using Test, SafeTestsets

@time begin
@time @safetestset "Third body force model tests..." begin include("thirdBodyPerterbationTests.jl") end
@time @safetestset "MEE tests..." begin include("meeTests.jl") end
end
