using Arblib
using TMatrix
using Test

const RTOL = 1e-6
const ATOL = 1e-6

@testset "TMatrix.jl" begin
    include("fixed.jl")
    include("fixed_arb.jl")
    include("random.jl")
end
