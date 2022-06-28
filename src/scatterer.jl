const R_1_2 = 1 // 2
const R_3_2 = 3 // 2
const R_1_4 = 1 // 4
const R_3_4 = 3 // 4
const R_1_3 = 1 // 3
const R_2_3 = 2 // 3
const R_4_3 = 4 // 3
const R_1_6 = 1 // 6

const DEFAULT_NCAP = Ref{Int64}(100)
const DEFAULT_NGCAP = Ref{Int64}(500)
const DEFAULT_NGCHEB = Ref{Int64}(60)
const ARB_APPROX_INV = Ref{Bool}(false)
const COLLECT_ACCURACY_INFO = Ref{Bool}(false)
const ACCURACY_INFO = Ref{Array{Tuple{String,Int64,Int64,Int64,Int64,Int64}}}([])
const FORCE_GC = Ref{Bool}(false)

const Q_INVERSION_WARNING = "Failed to invert Q numerically, try increasing precision or ngauss."

@doc raw"""
```
set_arb_approx_inv(t::Bool = true)
```

Turn on/off Arblib's approx_inv for ArbMatrix and AcbMatrix.
"""
function set_arb_approx_inv(t::Bool = true)
    ARB_APPROX_INV[] = t
    return
end

@doc raw"""
```
set_default_ncap(ncap::Int64 = 100)
```

Change the default `ncap` value.
"""
function set_default_ncap(ncap::Int64 = 100)
    DEFAULT_NCAP[] = ncap
    return
end

@doc raw"""
```
set_default_ngcap(ngcap::Int64 = 500)
```

Change the default `ngcap` value.
"""
function set_default_ngcap(ngcap::Int64 = 500)
    DEFAULT_NGCAP[] = ngcap
    return
end

@doc raw"""
```
set_default_ngcheb(ngcheb::Int64 = 60)
```

Change the default `ngcheb` value.
"""
function set_default_ngcheb(ngcheb::Int64 = 60)
    return DEFAULT_NGCHEB[] = ngcheb
end

function collect_accuracy_info(t::Bool = true)
    return COLLECT_ACCURACY_INFO[] = t
end

function clear_accuracy_info()
    TMatrix.ACCURACY_INFO[] = []
    return "Historical accuracy info has been cleared."
end

function save_accuracy_info()
    ts = format(now(UTC), "yyyy-mm-ddTHH:MM:SS")
    filename = "tmatrix_accuracy_info_" * ts * ".csv"
    return save_accuracy_info(filename)
end

function save_accuracy_info(filename::String)
    if COLLECT_ACCURACY_INFO[]
        open(filename, "w") do io
            print(io, "item,m,nmax,ngauss,original_accuracy,accuracy\n")
            for (item, m, nmax, ngauss, original_accuracy, accuracy) in ACCURACY_INFO[]
                print(io, item, ",", m, ",", nmax, ",", ngauss, ",", original_accuracy, ",", accuracy, "\n")
            end
        end
        return "Accuracy info saved to $filename"
    end

    return "Accuracy info is not collected."
end

@doc raw"""
Accompanied information of a scatterer.
"""
mutable struct ScattererInfo{RV<:AbstractVector,RM<:AbstractMatrix,CV<:AbstractVector,CM<:AbstractMatrix}
    nmax::Int64
    ngauss::Int64
    ncap::Int64
    ngcap::Int64
    an::RV # (nmax,)
    ann::RM # (nmax, nmax)
    sig::RV # (nmax,)
    x::RV # (ngauss,)
    w::RV # (ngauss,)
    s::RV # (ngauss,)
    r::RV # (ngauss,)
    dr::RV # (ngauss,)
    drr::RV # (ngauss,)
    wr²::RV # (ngauss,)
    kr⁻¹::RV # (ngauss,)
    kₛr⁻¹::CV # (ngauss,)
    d::RM # (ngauss, nmax)
    τ::RM # (ngauss, nmax)
    p::RM # (ngauss, nmax)
    jkr::RM # (ngauss, nmax)
    djkr::RM # (ngauss, nmax)
    j_tmp::RV # (2ncap,)
    ykr::RM # (ngauss, nmax)
    dykr::RM # (ngauss, nmax)
    hkr::CM # (ngauss, nmax)
    dhkr::CM # (ngauss, nmax)
    jkₛr::CM # (ngauss, nmax)
    djkₛr::CM # (ngauss, nmax)
    j_s_tmp::CV # (2ncap,)
    J₁₁::CM # (nmax, nmax)
    J₁₂::CM # (nmax, nmax)
    J₂₁::CM # (nmax, nmax)
    J₂₂::CM # (nmax, nmax)
    RgJ₁₁::CM # (nmax, nmax)
    RgJ₁₂::CM # (nmax, nmax)
    RgJ₂₁::CM # (nmax, nmax)
    RgJ₂₂::CM # (nmax, nmax)
    Q::CM # (2nmax, 2nmax)
    RgQ::CM # (2nmax, 2nmax)
end

@doc raw"""
Constructor of `ScattererInfo` for general data types. Space is pre-allocated to reduce allocations.
"""
function ScattererInfo(T::Type{<:Real})
    return ScattererInfo(
        0,
        0,
        DEFAULT_NCAP[],
        DEFAULT_NGCAP[],
        calc_an(T, DEFAULT_NCAP[]),
        calc_ann(T, DEFAULT_NCAP[]),
        calc_sig(T, DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[]),
        zeros(Complex{T}, DEFAULT_NGCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, 0),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(T, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NGCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, 0),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, DEFAULT_NCAP[], DEFAULT_NCAP[]),
        zeros(Complex{T}, 2DEFAULT_NCAP[], 2DEFAULT_NCAP[]),
        zeros(Complex{T}, 2DEFAULT_NCAP[], 2DEFAULT_NCAP[]),
    )
end

@doc raw"""
Abstract type for all scatterers.
"""
abstract type AbstractScatterer{T<:Real,CT<:Number} end

@doc raw"""
```
tmatrix_routine_mishchenko(scatterer::AbstractScatterer{T}, ddelta::T, ndgs::Int64) where {T<:Real}
```

Generate the routine function which follows the iteration strategy used by M. Mishchenko.

Parameters:

- `scatterer`: The scatterer.
- `ddelta`: The convergence threshold.
- `ndgs`: Determines the initial value of `ngauss` as `ndgs * nmax`. The initial value of `nmax` is determined by the formula $\max(4, \lfloor kr + 4.05 * \sqrt[3]{kr}\rfloor)$.
- `nstart`: Optional for manual setting of the initial value of `nmax`.
- `ngstart`: Optional for manual setting of the initial value of `ngauss`. This parameter takes effect only if `nstart≠0`.

"""
function tmatrix_routine_mishchenko(
    scatterer::AbstractScatterer{T},
    ddelta::T,
    ndgs::Int64;
    nstart::Int64 = 0,
    ngstart::Int64 = nstart * ndgs,
) where {T<:Real}
    kr = 2 * T(π) * volume_equivalent_radius(scatterer) / scatterer.λ
    if nstart == 0
        nstart = max(4, Int64(floor(kr + 4.05 * ∛kr)))
        ngstart = nstart * ndgs
    end
    Qsca0 = zero(T)
    Qext0 = zero(T)
    nmax_convergence = false

    function routine(ngauss::Int64, nmax::Int64, Qsca::T, Qext::T)
        ΔQsca = abs((Qsca0 - Qsca) / Qsca)
        ΔQext = abs((Qext0 - Qext) / Qext)
        Qsca0 = Qsca
        Qext0 = Qext
        Δ = max(ΔQsca, ΔQext)

        if !nmax_convergence
            @debug "[mishchenko] nmax iteration: nmax=$nmax ngauss=$ngauss ΔQsca=$ΔQsca ΔQext=$ΔQext"
            if Δ < ddelta
                @debug "[mishchenko] nmax convergence reached"
                nmax_convergence = true
                return ngauss + 2, nmax
            else
                return ngauss + ndgs, nmax + 1
            end
        else
            @debug "[mishchenko] ngauss iteration: nmax=$nmax ngauss=$ngauss ΔQsca=$ΔQsca ΔQext=$ΔQext"
            if Δ < ddelta
                @debug "[mishchenko] ngauss convergence reached"
                return -1, -1
            else
                return ngauss + 2, nmax
            end
        end
    end

    return ngstart, nstart, routine
end

@doc raw"""
```
tmatrix_routine_mishchenko_nmaxonly(scatterer::AbstractScatterer{T}, ddelta::T, ndgs::Int64) where {T<:Real}
```

Generate the routine function which generally follows the iteration strategy used by M. Mishchenko, but does not modify the value of `ngauss`.

Parameters:

- `scatterer`: The scatterer.
- `ddelta`: The convergence threshold.
- `ndgs`: Determines the initial value of `ngauss` as `ndgs * nmax`. The initial value of `nmax` is determined by the formula $\max(4, \lfloor kr + 4.05 * \sqrt[3]{kr}\rfloor)$.
- `nstart`: Optional for manual setting of the initial value of `nmax`.
- `ngstart`: Optional for manual setting of the initial value of `ngauss`. This parameter takes effect only if `nstart≠0`.

"""
function tmatrix_routine_mishchenko_nmaxonly(
    scatterer::AbstractScatterer{T},
    ddelta::T,
    ndgs::Int64;
    nstart::Int64 = 0,
    ngstart::Int64 = nstart * ndgs,
) where {T<:Real}
    kr = 2 * T(π) * volume_equivalent_radius(scatterer) / scatterer.λ
    if nstart == 0
        nstart = max(4, Int64(floor(kr + 4.05 * ∛kr)))
        ngstart = nstart * ndgs
    end
    Qsca0 = zero(T)
    Qext0 = zero(T)

    function routine(ngauss::Int64, nmax::Int64, Qsca::T, Qext::T)
        ΔQsca = abs((Qsca0 - Qsca) / Qsca)
        ΔQext = abs((Qext0 - Qext) / Qext)
        Qsca0 = Qsca
        Qext0 = Qext
        Δ = max(ΔQsca, ΔQext)

        @debug "[nmax-only mishchenko] nmax iteration: nmax=$nmax ngauss=$ngauss ΔQsca=$ΔQsca ΔQext=$ΔQext"
        if Δ < ddelta
            @debug "[nmax-only mishchenko] nmax convergence reached"
            return -1, -1
        else
            return ngauss, nmax + 1
        end
    end

    return ngstart, nstart, routine
end

@doc raw"""
```
calc_tmatrix!(scatterer::AbstractScatterer{T}) where {T<:Real}
```

Calculate the T-Matrix of the scatterer with default settings:

- Use `tmatrix_routine_mishchenko` to generate the routine function.
- Use `ddelta = 0.001` and `ndgs = 4`.

Parameters:

- `scatterer`: The scatterer.

"""
function calc_tmatrix!(scatterer::AbstractScatterer{T}) where {T<:Real}
    return calc_tmatrix!(scatterer, T(1) / T(1000), 4)
end

function calc_tmatrix!(scatterer::AbstractScatterer{T}, t::Tuple{Int64,Int64,Function};) where {T<:Real}
    ngauss, nmax, routine = t
    return calc_tmatrix!(scatterer, ngauss, nmax, routine)
end

function calc_tmatrix!(scatterer::AbstractScatterer{T}, ddelta::T, ndgs::Int64;) where {T<:Real}
    ngauss, nmax, routine = tmatrix_routine_mishchenko(scatterer, ddelta, ndgs)
    return calc_tmatrix!(scatterer, ngauss, nmax, routine)
end

@doc raw"""
```
calc_tmatrix(scatterer::AbstractScatterer{T}, ngstart::Int64, nstart::Int64, routine::Function) where {T<:Real}
```

Calculate the T-Matrix of the scatterer.

Parameters:

- `scatterer`: The scatterer.
- `ngstart`: The starting point of `ngauss`.
- `nstart`: The starting point of `nmax`.
- `routine`: The iteration routine function generated by a routine generator, internally or customly implemented.

"""
function calc_tmatrix!(
    scatterer::AbstractScatterer{T},
    ngstart::Int64,
    nstart::Int64,
    routine::Function,
) where {T<:Real}
    if T <: Arb
        @debug clear_accuracy_info()
    end

    ngauss, nmax = ngstart, nstart
    while true
        T0, _ = tmatr0!(scatterer, ngauss, nmax)
        Qext = sum((2n + 1) * real(T0[n, n] + T0[n + nmax, n + nmax]) for n in 1:nmax)
        Qsca = sum(
            (2n + 1) * real(T0[n, n] * T0[n, n]' + T0[n + nmax, n + nmax] * T0[n + nmax, n + nmax]') for n in 1:nmax
        )
        @debug "Qsca = $Qsca, Qext = $Qext"
        nngauss, nnmax = routine(ngauss, nmax, Qsca, Qext)
        if nnmax == -1
            return calc_tmatrix!(scatterer, ngauss, nmax, [T0])
        else
            ngauss, nmax = nngauss, nnmax
        end

        if FORCE_GC[]
            GC.gc() # Trigger GC manually to avoid memory leak
        end
    end
end

function calc_tmatrix!(
    scatterer::AbstractScatterer{T},
    ngauss::Int64,
    nmax::Int64,
    TT::Vector{<:AbstractMatrix},
) where {T<:Real}
    if length(TT) == 0
        @debug "Calculate T-Matrix for m = 0"
        T0, _ = tmatr0!(scatterer, ngauss, nmax)
        push!(TT, T0)
        if FORCE_GC[]
            GC.gc() # Trigger GC manually to avoid memory leak
        end
    end

    for m in 1:nmax
        @debug "Calculate T-Matrix for m = $m"
        Tm, _ = tmatr!(scatterer, m, ngauss, nmax)
        push!(TT, Tm)
        if FORCE_GC[]
            GC.gc() # Trigger GC manually to avoid memory leak
        end
    end

    @debug "Cross section" cross_section(TT, scatterer.λ)

    if T <: Arb
        @debug save_accuracy_info()
    end

    return TT
end

function calc_tmatrix!(scatterer::AbstractScatterer{T,CT}, ngauss::Int64, nmax::Int64) where {T<:Real,CT<:Number}
    @debug "Calculate T-Matrix for m = 0"
    T0, _ = tmatr0!(scatterer, ngauss, nmax)
    if FORCE_GC[]
        GC.gc() # Trigger GC manually to avoid memory leak
    end

    return calc_tmatrix!(scatterer, ngauss, nmax, [T0])
end

@doc raw"""
```
calc_amplitude(scatterer::AbstractScatterer{T}, α::T, β::T, ϑ_i::T, φ_i::T, ϑ_s::T, φ_s::T, TT::Vector{<:AbstractMatrix}) where {T<:Real}
```

Calculate the amplitude matrix and the phase matrix, given the scatterer and the geometry of the incident and the scattered beam. Use pre-computed T-Matrix when possible.

Parameters:

- `scatterer`: The scatterer.
- `α, β`: The Euler angle.
- `ϑ_i`: The zenith angle of the incident beam.
- `ϑ_s`: The zenith angle of the scattered beam.
- `φ_i`: The azimuth angle of the indicent beam.
- `φ_s`: The azimuth angle of the scatterer beam.
- `TT`: The pre-computed T-Matrix of the scatterer.

> All the angles here are input in degrees.
"""
function calc_amplitude(
    scatterer::AbstractScatterer{T},
    α::T,
    β::T,
    ϑ_i::T,
    ϑ_s::T,
    φ_i::T,
    φ_s::T,
    TT::Vector{<:AbstractMatrix},
) where {T<:Real}
    # Validate the input angles
    @assert 0.0 <= α <= 360.0 &&
            0.0 <= β <= 180.0 &&
            0.0 <= ϑ_i <= 180.0 &&
            0.0 <= ϑ_s <= 180.0 &&
            0.0 <= φ_i <= 360.0 &&
            0.0 <= φ_s <= 360.0

    nmax = length(TT) - 1

    α *= π / 180.0
    β *= π / 180.0
    ϑ_i *= π / 180.0
    ϑ_s *= π / 180.0
    φ_i *= π / 180.0
    φ_s *= π / 180.0

    cosβ = cos(β)
    sinβ = sin(β)

    cosϑ_i = cos(ϑ_i)
    sinϑ_i = sin(ϑ_i)
    cosφ = cos(φ_i - α)
    sinφ = sin(φ_i - α)
    cosϑ_p = cosϑ_i * cosβ + sinϑ_i * sinβ * cosφ
    ϑ_p = acos(cosϑ_p)
    cosφ_p = sinϑ_i * cosβ * cosφ - cosϑ_i * sinβ
    sinφ_p = sinϑ_i * sinφ
    φ_p = atan(sinφ_p, cosφ_p)

    cosϑ_s = cos(ϑ_s)
    sinϑ_s = sin(ϑ_s)
    cosφ = cos(φ_s - α)
    sinφ = sin(φ_s - α)
    cosϑ_q = cosϑ_s * cosβ + sinϑ_s * sinβ * cosφ
    ϑ_q = acos(cosϑ_q)
    cosφ_q = sinϑ_s * cosβ * cosφ - cosϑ_s * sinβ
    sinφ_q = sinϑ_s * sinφ
    φ_q = atan(sinφ_q, cosφ_q)

    B = zeros(T, 3, 3)
    cosα = cos(α)
    sinα = sin(α)
    B[1, 1] = cosα * cosβ
    B[1, 2] = sinα * cosβ
    B[1, 3] = -sinβ
    B[2, 1] = -sinα
    B[2, 2] = cosα
    B[2, 3] = 0.0
    B[3, 1] = cosα * sinβ
    B[3, 2] = sinα * sinβ
    B[3, 3] = cosβ

    AL = zeros(T, 3, 2)
    AL1 = zeros(T, 3, 2)
    cosφ_i = cos(φ_i)
    sinφ_i = sin(φ_i)
    cosφ_s = cos(φ_s)
    sinφ_s = sin(φ_s)
    AL[1, 1] = cosϑ_i * cosφ_i
    AL[1, 2] = -sinφ_i
    AL[2, 1] = cosϑ_i * sinφ_i
    AL[2, 2] = cosφ_i
    AL[3, 1] = -sinϑ_i
    AL1[1, 1] = cosϑ_s * cosφ_s
    AL1[1, 2] = -sinφ_s
    AL1[2, 1] = cosϑ_s * sinφ_s
    AL1[2, 2] = cosφ_s
    AL1[3, 1] = -sinϑ_s

    AP = zeros(T, 2, 3)
    AP1 = zeros(T, 2, 3)
    sinϑ_p = sin(ϑ_p)
    cosφ_p = cos(φ_p)
    sinφ_p = sin(φ_p)
    sinϑ_q = sin(ϑ_q)
    cosφ_q = cos(φ_q)
    sinφ_q = sin(φ_q)
    AP[1, 1] = cosϑ_p * cosφ_p
    AP[1, 2] = cosϑ_p * sinφ_p
    AP[1, 3] = -sinϑ_p
    AP[2, 1] = -sinφ_p
    AP[2, 2] = cosφ_p
    AP1[1, 1] = cosϑ_q * cosφ_q
    AP1[1, 2] = cosϑ_q * sinφ_q
    AP1[1, 3] = -sinϑ_q
    AP1[2, 1] = -sinφ_q
    AP1[2, 2] = cosφ_q

    R = AP * (B * AL)
    R1 = AP1 * (B * AL1)
    D = 1.0 / (R1[1, 1] * R1[2, 2] - R1[1, 2] * R1[2, 1])
    R1 = D * [R1[2, 2] -R1[1, 2]; -R1[2, 1] R1[1, 1]]

    CAL = [
        Complex{T}((1.0im)^(j - i - 1) * √((2j + 1) * (2i + 1) / (i * j * (i + 1) * (j + 1)))) for i in 1:nmax,
        j in 1:nmax
    ]

    φ = φ_q - φ_p
    VV = Complex{T}(0.0im)
    VH = Complex{T}(0.0im)
    HV = Complex{T}(0.0im)
    HH = Complex{T}(0.0im)
    for m in 0:nmax
        dv1, dv2 = vigampl(nmax, m, cosϑ_q)
        dv01, dv02 = vigampl(nmax, m, cosϑ_p)
        fc = 2.0cos(m * φ)
        fs = 2.0sin(m * φ)
        nm = nmax - max(m, 1) + 1
        dm = nmax - nm
        for nn in 1:nm
            for n in 1:nm
                T11 = TT[m + 1][n, nn]
                T22 = TT[m + 1][n + nm, nn + nm]
                if m == 0
                    CN = CAL[n + dm, nn + dm] * dv2[n + dm] * dv02[nn + dm]
                    VV += CN * T22
                    HH += CN * T11
                else
                    T12 = TT[m + 1][n, nn + nm]
                    T21 = TT[m + 1][n + nm, nn]
                    Cn₁ = CAL[n + dm, nn + dm] * fc
                    Cn₂ = CAL[n + dm, nn + dm] * fs
                    D11 = m^2 * dv1[n + dm] * dv01[nn + dm]
                    D12 = m * dv1[n + dm] * dv02[nn + dm]
                    D21 = m * dv2[n + dm] * dv01[nn + dm]
                    D22 = dv2[n + dm] * dv02[nn + dm]
                    VV += (T11 * D11 + T21 * D21 + T12 * D12 + T22 * D22) * Cn₁
                    VH += (T11 * D12 + T21 * D22 + T12 * D11 + T22 * D21) * Cn₂
                    HV -= (T11 * D21 + T21 * D11 + T12 * D22 + T22 * D12) * Cn₂
                    HH += (T11 * D22 + T21 * D12 + T12 * D21 + T22 * D11) * Cn₁
                end
            end
        end
    end

    DK = 2 * T(π) / scatterer.λ
    VV /= DK
    VH /= DK
    HV /= DK
    HH /= DK
    CVV = VV * R[1, 1] + VH * R[2, 1]
    CVH = VV * R[1, 2] + VH * R[2, 2]
    CHV = HV * R[1, 1] + HH * R[2, 1]
    CHH = HV * R[1, 2] + HH * R[2, 2]
    VV = R1[1, 1] * CVV + R1[1, 2] * CHV
    VH = R1[1, 1] * CVH + R1[1, 2] * CHH
    HV = R1[2, 1] * CVV + R1[2, 2] * CHV
    HH = R1[2, 1] * CVH + R1[2, 2] * CHH

    return [VV VH; HV HH]
end

calc_S = calc_amplitude

@doc raw"""
```
calc_phase(S::AbstractMatrix)
```

Calculate the phase matrix using the given amplitude matrix $\mathbf{S}$.

Parameters:

- `S`, the precalculated amplitude matrix.
"""
function calc_phase(S::AbstractMatrix)
    @assert size(S) == (2, 2)

    Z = zeros(eltype(S), 4, 4)
    Z[1, 1] = 0.5 * (S[1, 1] * S[1, 1]' + S[1, 2] * S[1, 2]' + S[2, 1] * S[2, 1]' + S[2, 2] * S[2, 2]')
    Z[1, 2] = 0.5 * (S[1, 1] * S[1, 1]' - S[1, 2] * S[1, 2]' + S[2, 1] * S[2, 1]' - S[2, 2] * S[2, 2]')
    Z[1, 3] = -S[1, 1] * S[1, 2]' - S[2, 2] * S[2, 1]'
    Z[1, 4] = 1.0im * (S[1, 1] * S[1, 2]' - S[2, 2] * S[2, 1]')
    Z[2, 1] = 0.5 * (S[1, 1] * S[1, 1]' + S[1, 2] * S[1, 2]' - S[2, 1] * S[2, 1]' - S[2, 2] * S[2, 2]')
    Z[2, 2] = 0.5 * (S[1, 1] * S[1, 1]' - S[1, 2] * S[1, 2]' - S[2, 1] * S[2, 1]' + S[2, 2] * S[2, 2]')
    Z[2, 3] = -S[1, 1] * S[1, 2]' + S[2, 2] * S[2, 1]'
    Z[2, 4] = 1.0im * (S[1, 1] * S[1, 2]' + S[2, 2] * S[2, 1]')
    Z[3, 1] = -S[1, 1] * S[2, 1]' - S[2, 2] * S[1, 2]'
    Z[3, 2] = -S[1, 1] * S[2, 1]' + S[2, 2] * S[1, 2]'
    Z[3, 3] = S[1, 1] * S[2, 2]' + S[1, 2] * S[2, 1]'
    Z[3, 4] = -1.0im * (S[1, 1] * S[2, 2]' + S[2, 1] * S[1, 2]')
    Z[4, 1] = 1.0im * (S[2, 1] * S[1, 1]' + S[2, 2] * S[1, 2]')
    Z[4, 2] = 1.0im * (S[2, 1] * S[1, 1]' - S[2, 2] * S[1, 2]')
    Z[4, 3] = -1.0im * (S[2, 2] * S[1, 1]' - S[1, 2] * S[2, 1]')
    Z[4, 4] = S[2, 2] * S[1, 1]' - S[1, 2] * S[2, 1]'

    return real.(Z)
end

@doc raw"""
```
calc_Z(S::AbstractMatrix)
```

Alias of `calc_phase`
"""
calc_Z = calc_phase

@doc raw"""
```
calc_SZ(scatterer::AbstractScatterer{T}, α::T, β::T, ϑ_i::T, φ_i::T, ϑ_s::T, φ_s::T, TT::Vector{<:AbstractMatrix}) where {T<:Real}
```

Calculate the S matrix and the Z matrix sequentially.
"""
function calc_SZ(
    scatterer::AbstractScatterer,
    α::T,
    β::T,
    ϑ_i::T,
    ϑ_s::T,
    φ_i::T,
    φ_s::T,
    TT::Vector{<:AbstractMatrix},
) where {T<:Real}
    S = calc_S(scatterer, α, β, ϑ_i, ϑ_s, φ_i, φ_s, TT)
    Z = calc_Z(S)

    return S, Z
end

@doc raw"""
```
hovenr(α₁, α₂, α₃, α₄, β₁, β₂)
```

Validate the expansion coefficients with the test of Van der Mee & Hovenier.

Note that the index of the coefficients should start from $0$.
"""
function hovenr(α₁, α₂, α₃, α₄, β₁, β₂)
    valid = true

    for l in eachindex(α₁)
        dl = 2l + 1
        ddl = 0.48dl
        if (l >= 1 && abs(α₁[l]) >= dl) ||
           abs(α₂[l]) >= dl ||
           abs(α₃[l]) >= dl ||
           abs(α₄[l]) >= dl ||
           abs(β₁[l]) >= ddl ||
           abs(β₂[l]) >= ddl
            @warn "Test of Van der Mee & Hovenier failed at L = $l"
            valid = false
            break
        end

        for c in 0:0.1:1.0
            cc = c^2
            c1 = cc * β₂[l]^2
            c2 = c * α₄[l]
            c3 = c * α₃[l]
            if (dl - c * α₁[l]) * (dl - c * α₂[l]) - cc * β₁[l]^2 <= -1e-4 ||
               (dl - c2) * (dl - c3) + c1 <= -1e-4 ||
               (dl + c2) * (dl - c3) - c1 <= -1e-4 ||
               (dl - c2) * (dl + c3) - c1 <= -1e-4
                valid = false
                break
            end
        end

        if !valid
            @warn "Test of Van der Mee & Hovenier failed at L = $l"
            break
        end
    end

    if valid
        @debug "Test of Van der Mee & Hovenier passed"
    end
end

@doc raw"""
```
calc_expansion_coefficients(TT::Vector{<:AbstractMatrix}, Csca::Real, λ::Real)
```

Calculate the expansion coefficients from a given T-Matrix.

Parameters:

- `TT`: The precalculated T-Matrix of a scatterer.
- `Csca`: The scattering cross setction.
- `λ`: The wavelength.

`Csca` and `λ` should have the same unit of length. E.g., if `λ` is in `μm`, `Csca` should be in `μm²`.
"""
function calc_expansion_coefficients(TT::Vector{<:AbstractMatrix}, Csca::Real, λ::Real)
    CT = eltype(TT[1])
    T = typeof(real(TT[1][1, 1]))
    nmax = length(TT) - 1
    Csca = T(Csca)
    λ = T(λ)

    ci = OffsetArray([CT(1.0im)^i for i in (-nmax):nmax], (-nmax):nmax)
    s = OffsetArray([T(2i + 1) for i in 0:(2nmax)], 0:(2nmax))
    ss = sqrt.(s)
    sig = OffsetArray([1 - 2 * (i % 2) for i in 0:(4nmax)], 0:(4nmax))

    T1 = OffsetArray(zeros(CT, 2nmax + 1, nmax), (-nmax):nmax, 1:nmax)
    T2 = OffsetArray(zeros(CT, 2nmax + 1, nmax), (-nmax):nmax, 1:nmax)
    A1 = zeros(CT, nmax)
    A2 = zeros(CT, nmax)
    B1 = OffsetArray(zeros(CT, 2nmax + 1, 2nmax + 1, nmax), 0:(2nmax), (-nmax):nmax, 1:nmax)
    B2 = OffsetArray(zeros(CT, 2nmax + 1, 2nmax + 1, nmax), 0:(2nmax), (-nmax):nmax, 1:nmax)

    for n in 1:nmax
        # Calculate T1 and T2
        for n′ in 1:nmax
            for m in 0:min(n, n′)
                nₒ = max(0, m - 1)
                nₗ = nmax - nₒ
                T11 = TT[m + 1][n - nₒ, n′ - nₒ]
                T12 = TT[m + 1][n - nₒ, n′ - nₒ + nₗ]
                T21 = TT[m + 1][n - nₒ + nₗ, n′ - nₒ]
                T22 = TT[m + 1][n - nₒ + nₗ, n′ - nₒ + nₗ]
                T1[m, n′] = T11 + T12 + T21 + T22
                T2[m, n′] = T11 + T12 - T21 - T22

                if m != 0
                    T1[-m, n′] = T11 - T12 - T21 + T22
                    T2[-m, n′] = T11 - T12 + T21 - T22
                end
            end
        end

        for n₁ in 0:(nmax + n)
            # Calculate A1 and A2
            for n′ in max(1, abs(n - n₁)):min(nmax, n₁ + n)
                A1[n′] = zero(CT)
                A2[n′] = zero(CT)
                for m₁ in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(T, n, m₁, n₁, 0, n′)
                    A1[n′] += cg * T1[m₁, n′]
                    A2[n′] += cg * T2[m₁, n′]
                end
                A1[n′] *= ci[n′ - n] / ss[n′]
                A2[n′] *= ci[n′ - n] / ss[n′]
            end

            # Calculate B1 and B2
            for m in max(1 - n₁, -n):min(n₁ + 1, n)
                for n′ in max(1, abs(n - n₁)):min(nmax, n₁ + n)
                    cg = clebschgordan(T, n, m, n₁, 1 - m, n′)
                    B1[n₁, m, n] += cg * A1[n′]
                    B2[n₁, m, n] += cg * A2[n′]
                end
            end
        end
    end

    # Calculate D
    D₀₀ = OffsetArray(zeros(T, 2nmax + 1, nmax, nmax), (-nmax):nmax, 1:nmax, 1:nmax)
    D₀₋₀ = OffsetArray(zeros(T, 2nmax + 1, nmax, nmax), (-nmax):nmax, 1:nmax, 1:nmax)
    D₂₂ = OffsetArray(zeros(T, 2nmax + 1, nmax, nmax), (-nmax):nmax, 1:nmax, 1:nmax)
    D₂₋₂ = OffsetArray(zeros(T, 2nmax + 1, nmax, nmax), (-nmax):nmax, 1:nmax, 1:nmax)
    D₀₂ = OffsetArray(zeros(CT, 2nmax + 1, nmax, nmax), (-nmax):nmax, 1:nmax, 1:nmax)

    for n in 1:nmax
        for n′ in 1:nmax
            for m in (-min(n, n′)):min(n, n′)
                for n₁ in abs(m - 1):(min(n, n′) + nmax)
                    D₀₀[m, n′, n] += s[n₁] * real(B1[n₁, m, n] * B1[n₁, m, n′]')
                    D₀₋₀[m, n′, n] += s[n₁] * real(B2[n₁, m, n] * B2[n₁, m, n′]')
                end
            end

            for m in max(-n, -n′ + 2):min(n, n′ + 2)
                for n₁ in abs(m - 1):(min(n, n′) + nmax)
                    D₂₂[m, n′, n] += s[n₁] * real(B1[n₁, m, n] * B1[n₁, 2 - m, n′]')
                    D₂₋₂[m, n′, n] += s[n₁] * real(B2[n₁, m, n] * B2[n₁, 2 - m, n′]')
                    D₀₂[m, n′, n] += s[n₁] * B2[n₁, m, n] * B1[n₁, 2 - m, n′]'
                end
            end
        end
    end

    h_const = λ^2 / (Csca * 4 * π)
    h = OffsetArray(
        [s[l] * h_const * ss[n] / ss[n′] for l in 0:(2nmax), n in 1:nmax, n′ in 1:nmax],
        0:(2nmax),
        1:nmax,
        1:nmax,
    )

    # Calculate g
    g₀₀ = OffsetArray(zeros(T, 2nmax + 1), 0:(2nmax))
    g₀₋₀ = OffsetArray(zeros(T, 2nmax + 1), 0:(2nmax))
    g₂₂ = OffsetArray(zeros(T, 2nmax + 1), 0:(2nmax))
    g₂₋₂ = OffsetArray(zeros(T, 2nmax + 1), 0:(2nmax))
    g₀₂ = OffsetArray(zeros(CT, 2nmax + 1), 0:(2nmax))

    for l in 0:(2nmax)
        for n in 1:nmax
            for n′ in max(1, abs(n - l)):min(nmax, n + l)
                cg1 = clebschgordan(T, n, 1, l, 0, n′)
                sm₀₀ = zero(T)
                sm₀₋₀ = zero(T)
                for m in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(T, n, m, l, 0, n′)
                    sm₀₀ += cg * D₀₀[m, n′, n]
                    sm₀₋₀ += cg * D₀₋₀[m, n′, n]
                end
                g₀₀[l] += h[l, n, n′] * cg1 * sm₀₀
                g₀₋₀[l] += h[l, n, n′] * cg1 * sig[n + n′ + l] * sm₀₋₀

                if l >= 2
                    cg2 = clebschgordan(T, n, -1, l, 2, n′)
                    sm₂₂ = zero(T)
                    sm₂₋₂ = zero(T)
                    sm₀₂ = zero(CT)
                    for m in max(-n, -n′ + 2):min(n, n′ + 2)
                        cg = clebschgordan(T, n, -m, l, 2, n′)
                        sm₂₂ += cg * D₂₂[m, n′, n]
                        sm₂₋₂ += cg * D₂₋₂[m, n′, n]
                        sm₀₂ += cg * D₀₂[m, n′, n]
                    end
                    g₂₂[l] += h[l, n, n′] * cg2 * sm₂₂
                    g₂₋₂[l] += h[l, n, n′] * cg2 * sig[n + n′ + l] * sm₂₋₂
                    g₀₂[l] += -h[l, n, n′] * cg1 * sm₀₂
                end
            end
        end
    end

    α₁ = g₀₀ + g₀₋₀
    α₂ = g₂₂ + g₂₋₂
    α₃ = g₂₂ - g₂₋₂
    α₄ = g₀₀ - g₀₋₀
    β₁ = 2real.(g₀₂)
    β₂ = 2imag.(g₀₂)

    # Validate the expansion coefficients
    @debug hovenr(α₁, α₂, α₃, α₄, β₁, β₂)

    return α₁, α₂, α₃, α₄, β₁, β₂
end

@doc raw"""
```
calc_scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θ)
```

Calculate the scatterering matrix elements from the given expansion coefficients.

Parameters:

- `α₁`, `α₂`, `α₃`, `α₄`, `β₁`, `β₂`: The precalculated expansion coefficients.
- `θ`: The scattering angle in degrees.
"""
function calc_scattering_matrix(
    α₁::AbstractVector{T},
    α₂::AbstractVector{T},
    α₃::AbstractVector{T},
    α₄::AbstractVector{T},
    β₁::AbstractVector{T},
    β₂::AbstractVector{T},
    θ::Real,
) where {T<:Real}
    lmax = length(α₁) - 1
    θ = Float64(θ) / 180 * π

    F₁₁ = sum(α₁[l] * WignerD.wignerdjmn(l, 0, 0, θ) for l in 0:lmax)
    F₂₂₊₃₃ = sum((α₂[l] + α₃[l]) * WignerD.wignerdjmn(l, 2, 2, θ) for l in 2:lmax)
    F₂₂₋₃₃ = sum((α₂[l] - α₃[l]) * WignerD.wignerdjmn(l, 2, -2, θ) for l in 2:lmax)
    F₂₂ = (F₂₂₊₃₃ + F₂₂₋₃₃) / 2
    F₃₃ = F₂₂₊₃₃ - F₂₂
    F₄₄ = sum(α₄[l] * WignerD.wignerdjmn(l, 0, 0, θ) for l in 0:lmax)
    F₁₂ = -sum(β₁[l] * WignerD.wignerdjmn(l, 0, 2, θ) for l in 2:lmax)
    F₃₄ = -sum(β₂[l] * WignerD.wignerdjmn(l, 0, 2, θ) for l in 2:lmax)

    return F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄
end

@doc raw"""
```
calc_scattering_matrix(scatterer, TT, Nθ)
```

Calculate the scatterering matrix elements from the given scatterer and precalculated T-Matrix.

Parameters:

- `scatterer`: The scatterer.
- `TT`: The T-Matrix.
- `Nθ`: Number of θ intervals (so the result will have `Nθ + 1` rows).
"""
function calc_scattering_matrix(scatterer::AbstractScatterer, TT::Vector{<:AbstractMatrix}, Nθ::Integer)
    Csca, _, _ = cross_section(TT, scatterer.λ)
    α₁, α₂, α₃, α₄, β₁, β₂ = calc_expansion_coefficients(TT, Csca, scatterer.λ)

    θ = Vector(range(0.0, 180.0, Nθ + 1))
    data = hcat(θ, reduce(hcat, map(collect, map(θᵢ -> calc_scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θᵢ), θ)))')
    df = DataFrame(data, ["θ", "s11", "s12", "s22", "s33", "s34", "s44"])
    return df
end

function theta_split(scatterer::AbstractScatterer{T}, ngauss::Int64) where {T<:Real}
    x = zeros(ngauss)
    w = zeros(ngauss)
    theta_split!(scatterer, ngauss, x, w)
    return x, w
end

@doc raw"""
```
calc_r(scatterer::Scatterer, ngauss::Int64)
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer.
"""
function calc_r(scatterer::AbstractScatterer{T}, ngauss::Int64) where {T<:Real}
    x = zeros(T, ngauss)
    w = zeros(T, ngauss)
    r = zeros(T, ngauss)
    dr = zeros(T, ngauss)
    calc_r!(scatterer, ngauss, x, w, r, dr)
    return r, dr
end

function update!(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64) where {T<:Real}
    info = scatterer.info

    # No need to recalculate if both `ngauss` and `nmax` remains the same.
    if ngauss == info.ngauss && nmax == info.nmax
        return
    end

    ncap = info.ncap
    ngcap = info.ngcap

    # Need to enlarge the arrays if either `ngauss` or `nmax` exceeds the current capacity.
    if ngauss > info.ngcap || nmax > info.ncap
        info.ncap = max(ncap, nmax * 2)
        info.ngcap = max(ngcap, ngauss * 2)

        if info.ncap > ncap
            info.an = calc_an(T, info.ncap)
            info.ann = calc_ann(T, info.ncap)
            info.sig = calc_sig(T, info.ncap)
        end

        if info.ngcap > ngcap
            info.x = zeros(T, info.ngcap)
            info.w = zeros(T, info.ngcap)
            info.s = zeros(T, info.ngcap)
            info.r = zeros(T, info.ngcap)
            info.dr = zeros(T, info.ngcap)
            info.drr = zeros(T, info.ngcap)
            info.wr² = zeros(T, info.ngcap)
            info.kr⁻¹ = zeros(T, info.ngcap)
            info.kₛr⁻¹ = zeros(Complex{T}, info.ngcap)
        end

        if info.ncap > ncap
            info.j_tmp = zeros(T, 3info.ncap)
            info.j_s_tmp = zeros(Complex{T}, 3info.ncap)
            info.J₁₁ = zeros(Complex{T}, info.ncap, info.ncap)
            info.J₁₂ = zeros(Complex{T}, info.ncap, info.ncap)
            info.J₂₁ = zeros(Complex{T}, info.ncap, info.ncap)
            info.J₂₂ = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ₁₁ = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ₁₂ = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ₂₁ = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ₂₂ = zeros(Complex{T}, info.ncap, info.ncap)
            info.Q = zeros(Complex{T}, 2info.ncap, 2info.ncap)
            info.RgQ = zeros(Complex{T}, 2info.ncap, 2info.ncap)
        end

        if info.ncap > ncap || info.ngcap > ngcap
            info.d = zeros(T, info.ngcap, info.ncap)
            info.τ = zeros(T, info.ngcap, info.ncap)
            info.p = zeros(T, info.ngcap, info.ncap)
            info.jkr = zeros(T, info.ngcap, info.ncap)
            info.djkr = zeros(T, info.ngcap, info.ncap)
            info.ykr = zeros(T, info.ngcap, info.ncap)
            info.dykr = zeros(T, info.ngcap, info.ncap)
            info.hkr = zeros(Complex{T}, info.ngcap, info.ncap)
            info.dhkr = zeros(Complex{T}, info.ngcap, info.ncap)
            info.jkₛr = zeros(Complex{T}, info.ngcap, info.ncap)
            info.djkₛr = zeros(Complex{T}, info.ngcap, info.ncap)
        end
    end

    # Need to recalculate `x`, `w`, `s`, `r`, `dr`, `kr⁻¹` and `kₛr⁻¹` only if `ngauss` changes.
    k = 2 * T(π) / scatterer.λ

    if ngauss != info.ngauss
        calc_r!(
            scatterer,
            ngauss,
            view(info.x, 1:ngauss),
            view(info.w, 1:ngauss),
            view(info.r, 1:ngauss),
            view(info.dr, 1:ngauss),
        )
    end

    kr = k * info.r[1:ngauss]
    kₛr = scatterer.m * kr

    if ngauss != info.ngauss
        @. info.drr[1:ngauss] = info.dr[1:ngauss] / info.r[1:ngauss]
        @. info.wr²[1:ngauss] = info.w[1:ngauss] * info.r[1:ngauss] * info.r[1:ngauss]
        @. info.s[1:ngauss] = 1 / (sin ∘ acos)(info.x[1:ngauss])
        @. info.kr⁻¹[1:ngauss] = 1 / kr
        @. info.kₛr⁻¹[1:ngauss] = 1 / kₛr
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.

    rmax = maximum(info.r[1:ngauss])
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(max(krmax, nmax) - nmax + 8.0 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb - nmax + 4.0 * ∛tb + 8.0 * √tb + 5))

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, nnmax1, view(info.jkr, i, :), view(info.djkr, i, :), info.j_tmp)
        sphericalbessely!(kr[i], nmax, view(info.ykr, i, :), view(info.dykr, i, :))
        sphericalbesselj!(kₛr[i], nmax, nnmax2, view(info.jkₛr, i, :), view(info.djkₛr, i, :), info.j_s_tmp)
    end

    view(info.hkr, 1:ngauss, 1:nmax) .= complex.(view(info.jkr, 1:ngauss, 1:nmax), view(info.ykr, 1:ngauss, 1:nmax))
    view(info.dhkr, 1:ngauss, 1:nmax) .= complex.(view(info.djkr, 1:ngauss, 1:nmax), view(info.dykr, 1:ngauss, 1:nmax))

    info.ngauss = ngauss
    info.nmax = nmax

    return
end

function constant(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64) where {T<:Real}
    an = [Float64(n * (n + 1)) for n in 1:nmax]
    ann = [0.5 * √((2n₁ + 1) * (2n₂ + 1) / (n₁ * (n₁ + 1) * n₂ * (n₂ + 1))) for n₁ in 1:nmax, n₂ in 1:nmax]

    x = zeros(ngauss)
    w = zeros(ngauss)
    theta_split!(scatterer, ngauss, x, w)
    s = 1 ./ (sin ∘ acos).(x)
    ss = s .^ 2

    return x, w, an, ann, s, ss
end

function vary(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64) where {T<:Real}
    r, dr = calc_r(scatterer, ngauss)
    λ = scatterer.λ
    k = 2 * T(π) / λ
    kr = k * r
    kr⁻¹ = 1.0 ./ kr
    kₛr = scatterer.m * kr
    kₛr⁻¹ = 1.0 ./ kₛr
    rmax = maximum(r)
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(8.0 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb + 4.0 * ∛tb + 8.0 * √tb - nmax + 5))

    jkr = zeros(T, ngauss, nmax)
    djkr = zeros(T, ngauss, nmax)
    ykr = zeros(T, ngauss, nmax)
    dykr = zeros(T, ngauss, nmax)
    jkₛr = zeros(Complex{T}, ngauss, nmax)
    djkₛr = zeros(Complex{T}, ngauss, nmax)
    j_tmp = zeros(T, nmax + nnmax1)
    j_s_tmp = zeros(Complex{T}, nmax + nnmax2)

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, nnmax1, view(jkr, i, :), view(djkr, i, :), j_tmp)
        sphericalbessely!(kr[i], nmax, view(ykr, i, :), view(dykr, i, :))
        sphericalbesselj!(kₛr[i], nmax, nnmax2, view(jkₛr, i, :), view(djkₛr, i, :), j_s_tmp)
    end

    return r, dr, kr⁻¹, kₛr⁻¹, jkr, djkr, ykr, dykr, jkₛr, djkₛr
end

function tmatr0!(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64;) where {T<:Real}
    @assert iseven(ngauss)

    CT = Complex{T}

    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    an = view(info.an, 1:nmax)
    ann = view(info.ann, 1:nmax, 1:nmax)
    sig = view(info.sig, 1:nmax)
    x = view(info.x, 1:ngauss)

    d = view(info.d, 1:ngauss, 1:nmax)
    τ = view(info.τ, 1:ngauss, 1:nmax)
    for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, 0, x[i], view(d, i, :), view(τ, i, :))
        d[ineg, :] .= view(d, i, :) .* sig
        τ[ineg, :] .= view(τ, i, :) .* sig .* (-1)
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    drr = view(info.drr, 1:ngss)
    kr⁻¹ = view(info.kr⁻¹, 1:ngss)
    kₛr⁻¹ = view(info.kₛr⁻¹, 1:ngss)
    wr² = view(info.wr², 1:ngss)

    jkr = view(info.jkr, 1:ngss, 1:nmax)
    djkr = view(info.djkr, 1:ngss, 1:nmax)
    hkr = view(info.hkr, 1:ngss, 1:nmax)
    dhkr = view(info.dhkr, 1:ngss, 1:nmax)
    jkₛr = view(info.jkₛr, 1:ngss, 1:nmax)
    djkₛr = view(info.djkₛr, 1:ngss, 1:nmax)

    J₁₂ = view(info.J₁₂, 1:nmax, 1:nmax)
    J₂₁ = view(info.J₂₁, 1:nmax, 1:nmax)
    RgJ₁₂ = view(info.RgJ₁₂, 1:nmax, 1:nmax)
    RgJ₂₁ = view(info.RgJ₂₁, 1:nmax, 1:nmax)
    fill!(J₁₂, zero(CT))
    fill!(J₂₁, zero(CT))
    fill!(RgJ₁₂, zero(CT))
    fill!(RgJ₂₁, zero(CT))

    Threads.@threads for nn in 0:(nmax * nmax - 1)
        n₂ = nn ÷ nmax + 1
        n₁ = nn % nmax + 1
        if !(sym && (n₁ + n₂) % 2 == 1)
            for i in 1:ngss
                τ₁τ₂ = τ[i, n₁] * τ[i, n₂]
                d₁τ₂ = d[i, n₁] * τ[i, n₂]
                d₂τ₁ = d[i, n₂] * τ[i, n₁]

                J₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (dhkr[i, n₁] * τ₁τ₂ + drr[i] * an[n₁] * hkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                J₂₁[n₁, n₂] +=
                    wr²[i] * hkr[i, n₁] * (djkₛr[i, n₂] * τ₁τ₂ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)

                RgJ₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (djkr[i, n₁] * τ₁τ₂ + drr[i] * an[n₁] * jkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                RgJ₂₁[n₁, n₂] +=
                    wr²[i] * jkr[i, n₁] * (djkₛr[i, n₂] * τ₁τ₂ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)
            end
        end
    end

    @. J₁₂ *= -1.0im * ann
    @. J₂₁ *= 1.0im * ann

    @. RgJ₁₂ *= -1.0im * ann
    @. RgJ₂₁ *= 1.0im * ann

    k = 2 * T(π) / scatterer.λ
    kₛ = k * scatterer.m
    k² = k^2
    kkₛ = k * kₛ

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    Q = view(info.Q, 1:(2nmax), 1:(2nmax))
    Q₁₁ = view(Q, 1:nmax, 1:nmax)
    Q₂₂ = view(Q, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    RgQ = view(info.RgQ, 1:(2nmax), 1:(2nmax))
    RgQ₁₁ = view(RgQ, 1:nmax, 1:nmax)
    RgQ₂₂ = view(RgQ, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    fill!(Q, zero(CT))
    fill!(RgQ, zero(CT))

    @. Q₁₁ = kkₛ * J₂₁ + k² * J₁₂
    @. Q₂₂ = kkₛ * J₁₂ + k² * J₂₁

    @. RgQ₁₁ = kkₛ * RgJ₂₁ + k² * RgJ₁₂
    @. RgQ₂₂ = kkₛ * RgJ₁₂ + k² * RgJ₂₁

    T0 = RgQ * inv(Q)
    T0 .*= -1

    return T0, Q, RgQ
end

function tmatr!(scatterer::AbstractScatterer{T}, m::Int64, ngauss::Int64, nmax::Int64;) where {T<:Real}
    @assert iseven(ngauss)

    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    mm = max(m, 1)
    an = view(info.an, 1:nmax)
    ann = view(info.ann, 1:nmax, 1:nmax)
    sig = view(info.sig, 1:nmax)
    x = view(info.x, 1:ngauss)
    s = view(info.s, 1:ngauss)

    d = view(info.d, 1:ngauss, 1:nmax)
    p = view(info.p, 1:ngauss, 1:nmax)
    τ = view(info.τ, 1:ngauss, 1:nmax)

    for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, m, x[i], view(d, i, :), view(τ, i, :))
        p[i, :] .= view(d, i, :) .* (s[i] * m)
        d[ineg, :] .= view(d, i, :) .* sig
        τ[ineg, :] .= view(τ, i, :) .* sig .* (-1.0)
        p[ineg, :] .= view(p, i, :) .* sig
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    drr = view(info.drr, 1:ngss)
    kr⁻¹ = view(info.kr⁻¹, 1:ngss)
    kₛr⁻¹ = view(info.kₛr⁻¹, 1:ngss)
    wr² = view(info.wr², 1:ngss)

    jkr = view(info.jkr, 1:ngss, 1:nmax)
    djkr = view(info.djkr, 1:ngss, 1:nmax)
    hkr = view(info.hkr, 1:ngss, 1:nmax)
    dhkr = view(info.dhkr, 1:ngss, 1:nmax)
    jkₛr = view(info.jkₛr, 1:ngss, 1:nmax)
    djkₛr = view(info.djkₛr, 1:ngss, 1:nmax)

    J₁₁ = view(info.J₁₁, mm:nmax, mm:nmax)
    J₁₂ = view(info.J₁₂, mm:nmax, mm:nmax)
    J₂₁ = view(info.J₂₁, mm:nmax, mm:nmax)
    J₂₂ = view(info.J₂₂, mm:nmax, mm:nmax)
    RgJ₁₁ = view(info.RgJ₁₁, mm:nmax, mm:nmax)
    RgJ₁₂ = view(info.RgJ₁₂, mm:nmax, mm:nmax)
    RgJ₂₁ = view(info.RgJ₂₁, mm:nmax, mm:nmax)
    RgJ₂₂ = view(info.RgJ₂₂, mm:nmax, mm:nmax)
    fill!(J₁₁, zero(eltype(J₁₁)))
    fill!(J₁₂, zero(eltype(J₁₂)))
    fill!(J₂₁, zero(eltype(J₂₁)))
    fill!(J₂₂, zero(eltype(J₂₂)))
    fill!(RgJ₁₁, zero(eltype(RgJ₁₁)))
    fill!(RgJ₁₂, zero(eltype(RgJ₁₂)))
    fill!(RgJ₂₁, zero(eltype(RgJ₂₁)))
    fill!(RgJ₂₂, zero(eltype(RgJ₂₂)))

    OffsetJ₁₁ = OffsetArray(J₁₁, mm:nmax, mm:nmax)
    OffsetJ₁₂ = OffsetArray(J₁₂, mm:nmax, mm:nmax)
    OffsetJ₂₁ = OffsetArray(J₂₁, mm:nmax, mm:nmax)
    OffsetJ₂₂ = OffsetArray(J₂₂, mm:nmax, mm:nmax)
    OffsetRgJ₁₁ = OffsetArray(RgJ₁₁, mm:nmax, mm:nmax)
    OffsetRgJ₁₂ = OffsetArray(RgJ₁₂, mm:nmax, mm:nmax)
    OffsetRgJ₂₁ = OffsetArray(RgJ₂₁, mm:nmax, mm:nmax)
    OffsetRgJ₂₂ = OffsetArray(RgJ₂₂, mm:nmax, mm:nmax)

    nm = nmax - mm + 1

    Threads.@threads for nn in 0:(nm * nm - 1)
        n₂ = nn ÷ nm + mm
        n₁ = nn % nm + mm
        if !(sym && (n₁ + n₂) % 2 == 0)
            for i in 1:ngss
                pττp = p[i, n₁] * τ[i, n₂] + p[i, n₂] * τ[i, n₁]
                p₁d₂ = p[i, n₁] * d[i, n₂]

                OffsetJ₁₁[n₁, n₂] += wr²[i] * hkr[i, n₁] * jkₛr[i, n₂] * pττp

                OffsetJ₂₂[n₁, n₂] +=
                    wr²[i] * (
                        dhkr[i, n₁] * djkₛr[i, n₂] * pττp +
                        drr[i] *
                        (an[n₁] * hkr[i, n₁] * kr⁻¹[i] * djkₛr[i, n₂] + an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * dhkr[i, n₁]) *
                        p₁d₂
                    )

                OffsetRgJ₁₁[n₁, n₂] += wr²[i] * jkr[i, n₁] * jkₛr[i, n₂] * pττp

                OffsetRgJ₂₂[n₁, n₂] +=
                    wr²[i] * (
                        djkr[i, n₁] * djkₛr[i, n₂] * pττp +
                        drr[i] *
                        (an[n₁] * jkr[i, n₁] * kr⁻¹[i] * djkₛr[i, n₂] + an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * djkr[i, n₁]) *
                        p₁d₂
                    )
            end
        end

        if !(sym && (n₁ + n₂) % 2 == 1)
            for i in 1:ngss
                ppττ = p[i, n₁] * p[i, n₂] + τ[i, n₁] * τ[i, n₂]
                d₁τ₂ = d[i, n₁] * τ[i, n₂]
                d₂τ₁ = d[i, n₂] * τ[i, n₁]

                OffsetJ₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (dhkr[i, n₁] * ppττ + drr[i] * an[n₁] * hkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                OffsetJ₂₁[n₁, n₂] +=
                    wr²[i] * hkr[i, n₁] * (djkₛr[i, n₂] * ppττ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)

                OffsetRgJ₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (djkr[i, n₁] * ppττ + drr[i] * an[n₁] * jkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                OffsetRgJ₂₁[n₁, n₂] +=
                    wr²[i] * jkr[i, n₁] * (djkₛr[i, n₂] * ppττ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)
            end
        end
    end

    ann = view(ann, mm:nmax, mm:nmax)
    @. J₁₁ *= -1 * ann
    @. J₁₂ *= -1im * ann
    @. J₂₁ *= 1im * ann
    @. J₂₂ *= -1 * ann

    @. RgJ₁₁ *= -1 * ann
    @. RgJ₁₂ *= -1im * ann
    @. RgJ₂₁ *= 1im * ann
    @. RgJ₂₂ *= -1 * ann

    k = 2 * T(π) / scatterer.λ
    kₛ = k * scatterer.m
    k² = k^2
    kkₛ = k * kₛ

    Q = view(info.Q, 1:(2nm), 1:(2nm))
    Q₁₁ = view(Q, 1:nm, 1:nm)
    Q₁₂ = view(Q, 1:nm, (nm + 1):(2nm))
    Q₂₁ = view(Q, (nm + 1):(2nm), 1:nm)
    Q₂₂ = view(Q, (nm + 1):(2nm), (nm + 1):(2nm))
    RgQ = view(info.RgQ, 1:(2nm), 1:(2nm))
    RgQ₁₁ = view(RgQ, 1:nm, 1:nm)
    RgQ₁₂ = view(RgQ, 1:nm, (nm + 1):(2nm))
    RgQ₂₁ = view(RgQ, (nm + 1):(2nm), 1:nm)
    RgQ₂₂ = view(RgQ, (nm + 1):(2nm), (nm + 1):(2nm))
    fill!(Q, zero(eltype(Q)))
    fill!(RgQ, zero(eltype(RgQ)))

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    @. Q₁₁ = kkₛ * J₂₁ + k² * J₁₂
    @. Q₁₂ = kkₛ * J₁₁ + k² * J₂₂
    @. Q₂₁ = kkₛ * J₂₂ + k² * J₁₁
    @. Q₂₂ = kkₛ * J₁₂ + k² * J₂₁

    @. RgQ₁₁ = kkₛ * RgJ₂₁ + k² * RgJ₁₂
    @. RgQ₁₂ = kkₛ * RgJ₁₁ + k² * RgJ₂₂
    @. RgQ₂₁ = kkₛ * RgJ₂₂ + k² * RgJ₁₁
    @. RgQ₂₂ = kkₛ * RgJ₁₂ + k² * RgJ₂₁

    Tm = RgQ * inv(Q)
    Tm .*= -1

    return Tm, Q, RgQ
end

function check_convergece(
    scatterer::AbstractScatterer{T},
    nmax::Int64,
    ngauss::Int64,
    ndgs::Int64,
    ddelta::Real,
) where {T<:Real}
    T0, _ = tmatr0!(scatterer, ngauss, nmax)
    Qext = sum((2n + 1) * real(T0[n, n] + T0[n + nmax, n + nmax]) for n in 1:nmax)
    Qsca = sum((2n + 1) * real(T0[n, n] * T0[n, n]' + T0[n + nmax, n + nmax] * T0[n + nmax, n + nmax]') for n in 1:nmax)

    T0′, _ = tmatr0!(scatterer, ngauss + ndgs, nmax + 1)
    Qext′ = sum((2n + 1) * real(T0′[n, n] + T0′[n + nmax, n + nmax]) for n in 1:nmax)
    Qsca′ = sum(
        (2n + 1) * real(T0′[n, n] * T0′[n, n]' + T0′[n + nmax, n + nmax] * T0′[n + nmax, n + nmax]') for n in 1:nmax
    )

    ΔQext = abs((Qext - Qext′) / abs(Qext′))
    ΔQsca = abs((Qsca - Qsca′) / abs(Qsca′))

    return max(ΔQext, ΔQsca) <= ddelta
end
