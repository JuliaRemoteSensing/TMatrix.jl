module Random

using OffsetArrays
using TMatrix_jll

const NPN1 = 100
const NPN2 = 2NPN1
const NPN3 = NPN1 + 1
const NPN4 = 80
const NPN5 = 2NPN4
const NPN6 = NPN4 + 1
const NPL = NPN2 + 1
const NPL1 = NPN5 + 1
const NPNG1 = 300
const NPNG2 = 2NPNG1

function __init__()
    ccall((:fact_, tmatrix_random_orientation), Cvoid, ())
    ccall((:signum_, tmatrix_random_orientation), Cvoid, ())
    return
end

function set_T!(T::Vector{<:AbstractMatrix})
    tmat_ptr = cglobal((:tmat_, tmatrix_random_orientation), UInt64)
    blk_len = NPN6 * NPN4 * NPN4 * 4
    TR11 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr), (NPN6, NPN4, NPN4))
    TR12 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len), (NPN6, NPN4, NPN4))
    TR21 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 2), (NPN6, NPN4, NPN4))
    TR22 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 3), (NPN6, NPN4, NPN4))
    TI11 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 4), (NPN6, NPN4, NPN4))
    TI12 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 5), (NPN6, NPN4, NPN4))
    TI21 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 6), (NPN6, NPN4, NPN4))
    TI22 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 7), (NPN6, NPN4, NPN4))

    fill!(TR11, 0.0f0)
    fill!(TR12, 0.0f0)
    fill!(TR21, 0.0f0)
    fill!(TR22, 0.0f0)
    fill!(TI11, 0.0f0)
    fill!(TI12, 0.0f0)
    fill!(TI21, 0.0f0)
    fill!(TI22, 0.0f0)

    nmax = length(T) - 1
    for m in 0:nmax
        m1 = nmax - max(0, m - 1)
        TR11[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ real).(T[m + 1][1:m1, 1:m1])
        TR12[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ real).(T[m + 1][1:m1, (m1 + 1):(2m1)])
        TR21[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ real).(T[m + 1][(m1 + 1):(2m1), 1:m1])
        TR22[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] =
            (Float32 ∘ real).(T[m + 1][(m1 + 1):(2m1), (m1 + 1):(2m1)])
        TI11[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ imag).(T[m + 1][1:m1, 1:m1])
        TI12[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ imag).(T[m + 1][1:m1, (m1 + 1):(2m1)])
        TI21[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] = (Float32 ∘ imag).(T[m + 1][(m1 + 1):(2m1), 1:m1])
        TI22[m + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax] =
            (Float32 ∘ imag).(T[m + 1][(m1 + 1):(2m1), (m1 + 1):(2m1)])
    end
end

function ccg(n::Int64, n1::Int64, nmax::Int64, k1::Int64, k2::Int64)
    gg = zeros(NPL1, NPN6)

    ccall(
        (:ccg_, tmatrix_random_orientation),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}),
        convert(Int32, n),
        convert(Int32, n1),
        convert(Int32, nmax),
        convert(Int32, k1),
        convert(Int32, k2),
        gg,
    )

    return gg[(NPN6 - n):(NPN6 + n), 1:(min(n + n1, nmax) + 1)]
end

function gsp(nmax::Int64, Csca::Real, λ::Real)
    Csca = Float64(Csca)
    λ = Float64(λ)

    α₁ = zeros(NPL)
    α₂ = zeros(NPL)
    α₃ = zeros(NPL)
    α₄ = zeros(NPL)
    β₁ = zeros(NPL)
    β₂ = zeros(NPL)
    lmax = Ref{Int32}(0)

    ccall(
        (:gsp_, tmatrix_random_orientation),
        Cvoid,
        (
            Ref{Int32},
            Ref{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Int32},
        ),
        convert(Int32, nmax),
        Csca,
        λ,
        α₁,
        α₂,
        α₃,
        α₄,
        β₁,
        β₂,
        lmax,
    )

    lmax = convert(Int64, lmax[])

    return OffsetArray(α₁[1:lmax], 0:(lmax - 1)),
    OffsetArray(α₂[1:lmax], 0:(lmax - 1)),
    OffsetArray(α₃[1:lmax], 0:(lmax - 1)),
    OffsetArray(α₄[1:lmax], 0:(lmax - 1)),
    OffsetArray(β₁[1:lmax], 0:(lmax - 1)),
    OffsetArray(β₂[1:lmax], 0:(lmax - 1)),
    lmax - 1
end

end
