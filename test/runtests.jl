ENV["NEMO_THREADED"] = 1

using Arblib
using TMatrix
using Test

Arblib.flint_set_num_threads(Threads.nthreads())
const RTOL = 1e-6
const ATOL = 1e-6

@testset "TMatrix.jl" begin
    include("fixed.jl")
    include("fixed_arb.jl")
    include("extra_shapes.jl")
    include("random.jl")
end
