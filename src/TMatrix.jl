module TMatrix

using Arblib
using DataFrames
using Dates: format, now, UTC
using FastGaussQuadrature
using LinearAlgebra
using OffsetArrays
using WignerD
using WignerSymbols

include("arb_compat.jl")
include("helpers.jl")
include("scatterer.jl")
include("scatterer_arb.jl")
include("shapes/shapes.jl")
include("wrapper/wrapper.jl")

export AbstractTMatrix, AxisymmetricTMatrix

abstract type AbstractTMatrix{T} <: AbstractArray{T,6} end

volume_equivalent_radius(::AbstractTMatrix) = error("Not implemented")

struct AxisymmetricTMatrix{T} <: AbstractTMatrix{T}
    internal::Vector{Matrix{T}}
    lmax::Int
    rev::T
    axes::NTuple{6,UnitRange{Int}}

    function AxisymmetricTMatrix(scatterer::AbstractScatterer, tm::Vector{Matrix{T}}) where {T}
        lmax = length(tm) - 1
        return new{T}(
            tm,
            lmax,
            volume_equivalent_radius(scatterer),
            (1:lmax, (-lmax):lmax, 1:2, 1:lmax, (-lmax):lmax, 1:2),
        )
    end
end

volume_equivalent_radius(tm::AxisymmetricTMatrix) = tm.rev
Base.size(tm::AxisymmetricTMatrix) = (tm.lmax, 2tm.lmax + 1, 2, tm.lmax, 2tm.lmax + 1, 2)
Base.axes(tm::AxisymmetricTMatrix) = tm.axes

function Base.getindex(tm::AxisymmetricTMatrix{T}, n::Int, m::Int, p::Int, n′::Int, m′::Int, q::Int) where {T}
    if m == m′ && abs(m) <= min(n, n′)
        ma = abs(m)
        n₀ = tm.lmax - max(1, ma) + 1
        nn = (3 - p) * n₀ + n - tm.lmax
        nn′ = (3 - q) * n₀ + n′ - tm.lmax
        tm.internal[ma + 1][nn, nn′]
    else
        zero(T)
    end
end

end
