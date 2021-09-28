@enum Shape SHAPE_SPHEROID SHAPE_CYLINDER SHAPE_CHEBYSHEV
@enum RadiusType RADIUS_EQUAL_VOLUME RADIUS_EQUAL_AREA RADIUS_MAXIMUM

const DEFAULT_NCAP = Ref{Int64}(100)
const DEFAULT_NGCAP = Ref{Int64}(500)
const DEFAULT_NGCHEB = Ref{Int64}(60)
const ARB_APPROX_INV = Ref{Bool}(false)
const COLLECT_ACCURACY_INFO = Ref{Bool}(false)
const ACCURACY_INFO = Ref{Array{Tuple{String,Int64,Int64,Int64,Int64,Int64}}}([])

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
    kr1::RV # (ngauss,)
    kr_s1::CV # (ngauss,)
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
    jkr_s::CM # (ngauss, nmax)
    djkr_s::CM # (ngauss, nmax)
    j_s_tmp::CV # (2ncap,)
    J11::CM # (nmax, nmax)
    J12::CM # (nmax, nmax)
    J21::CM # (nmax, nmax)
    J22::CM # (nmax, nmax)
    RgJ11::CM # (nmax, nmax)
    RgJ12::CM # (nmax, nmax)
    RgJ21::CM # (nmax, nmax)
    RgJ22::CM # (nmax, nmax)
    Q::CM # (2nmax, 2nmax)
    RgQ::CM # (2nmax, 2nmax)
end

@doc raw"""
Constructor of `ScattererInfo`.

For most types, space is pre-assigned to reduce allocations.

However, for `Arb`, pre-assignment is not used, since `SubArray` does not work harmoniously with `ArbMatrix`.
"""
function ScattererInfo(T)
    return T <: Arb ?
           ScattererInfo(
        0,
        0,
        0,
        0,
        ArbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        AcbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefVector(0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
    ) :
           ScattererInfo(
        0,
        0,
        DEFAULT_NCAP[],
        DEFAULT_NGCAP[],
        [T(n * (n + 1)) for n in 1:DEFAULT_NCAP[]],
        [
            T(0.5 * √((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1)))) for n1 in 1:DEFAULT_NCAP[],
            n2 in 1:DEFAULT_NCAP[]
        ],
        [i % 2 == 1 ? T(-1.0) : T(1.0) for i in 1:DEFAULT_NCAP[]],
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
A spheroid scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `a_to_c`: The ratio $a/c$ of the horizontal to rotational axes.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Spheroid{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    a_to_c::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

@doc raw"""
A cylinder scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `d_to_h`: The diameter-to-height ratio $D/H$.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Cylinder{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    d_to_h::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

@doc raw"""
A Chebyshev scatterer defined by

    $r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))$

in which $T_n(\cos\theta)=\cos n\theta$.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `ε`: The deformation parameter, which satisfies $0\le\varepsilon<1$.
- `n`: The degree of the Chebyshev polynomial.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Chebyshev{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    ε::T
    n::Int64
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

@doc raw"""
Scatterer constructor with named parameters.

Parameters:

- `T`: The data type to be used in future calculations.
- `r`: The equivalent radius.
- `shape`: The particle shape. Possible values are `SHAPE_SPHEROID` (default), `SHAPE_CYLINDER` and `SHAPE_CHEBYSHEV`.
- `axis_ratio`: For spheroids, it is the ratio $a/b$ of the horizontal to rotational axes. For cylinders, it is the diameter-to-length ratio $D/L$. For Chebyshev particles, it is the deformation parameter $\varepsilon$.
- `radius_type`: The type of the equivalent radius. Possible values are `RADIUS_EQUAL_VOLUME` (default), `RADIUS_EQUAL_AREA` and `RADIUS_MAXIMUM`. All radius types will be transformed into the equivalent volume radius.
- `refractive_index`: The complex refractive index.
- `ngauss`: Number of points for Gaussian integration. Required only for Chebyshev particles.
- `n`: Degree of the Chebyshev polynomial. Required only for Chebyshev particles.
- `λ`: The wavelength to be used in future calculations.

"""
function Scatterer(
    T::Type{<:Real} = Float64;
    r::Real = 1.0,
    shape::Shape = SHAPE_SPHEROID,
    axis_ratio::Real = 1.0,
    radius_type::RadiusType = RADIUS_EQUAL_VOLUME,
    refractive_index::Number = 1.0 + 0.0im,
    n::Int64 = 2,
    ngauss::Int64 = DEFAULT_NGCHEB[],
    λ::Real = 1.0,
)
    if shape == SHAPE_CHEBYSHEV && (abs(axis_ratio) < 0.0 || abs(axis_ratio) >= 1.0)
        error("Constraint violated: Chebyshev particles should have 0≤|ε|<1.")
    end

    r = T(r)
    axis_ratio = T(axis_ratio)
    refractive_index = Complex{T}(refractive_index)
    λ = T(λ)
    f12 = T(1 // 2)
    f32 = T(3 // 2)
    f14 = T(1 // 4)
    f34 = T(3 // 4)
    f13 = T(1 // 3)
    f23 = T(2 // 3)
    f43 = T(4 // 3)

    if radius_type == RADIUS_EQUAL_VOLUME
        rev = r
    elseif radius_type == RADIUS_EQUAL_AREA
        if shape == SHAPE_SPHEROID
            d = axis_ratio
            if abs(d - 1.0) < eps(d)
                ratio = 1.0
            elseif d > 1.0
                e = √(1 - 1 / d^2)
                ratio = 1 / √(f14 * (2d^f23 + d^(-f43) * log((1 + e) / (1 - e)) / e))
            elseif d < 1.0
                e = √(1 - d^2)
                ratio = 1 / √(f12 * (d^f23 + d^(-f13) * asin(e) / e))
            end
        elseif shape == SHAPE_CYLINDER
            e = axis_ratio
            ratio = ∛(f32 / e) / √((e + 2) / 2e)
        elseif shape == SHAPE_CHEBYSHEV
            e = axis_ratio
            en = e * n
            x, w = gausslegendre(T, ngauss)
            s = zero(T)
            v = zero(T)
            @simd for i in 1:ngauss
                θ = acos(x[i])
                nθ = n * θ
                sinθ = sin(θ)
                sinnθ = sin(nθ)
                cosnθ = cos(nθ)
                a = 1 + e * cosnθ
                ens = en * sinnθ
                s += w[i] * a * √(a^2 + ens^2)
                v += w[i] * (sinθ * a + x[i] * ens) * sinθ * a^2
            end
            rs = √(f12 * s)
            rv = ∛(f34 * v)
            ratio = rv / rs
        end

        rev = r * ratio
    elseif radius_type == RADIUS_MAXIMUM
        if shape == SHAPE_SPHEROID
            rev = axis_ratio > 1.0 ? r / ∛axis_ratio : r * axis_ratio^f23
        elseif shape == SHAPE_CYLINDER
            rev = axis_ratio > 1.0 ? r * ∛(f32 / axis_ratio) : r * ∛(f32 * axis_ratio^2)
        elseif shape == SHAPE_CHEBYSHEV
            e = axis_ratio
            en = e * n
            x, w = gausslegendre(T, ngauss)
            θ = acos.(x)
            nθ = n * θ
            sinθ = sin.(θ)
            sinnθ = sin.(nθ)
            cosnθ = cos.(nθ)
            a = e * cosnθ .+ 1
            ens = en * sinnθ
            v = sum(w .* (sinθ .* a + x .* ens) .* sinθ * a .^ 2)
            rv = ∛(f34 * v)
            rev = r / (1 + e) * rv
        end
    end

    rev = T(rev)

    if shape == SHAPE_SPHEROID
        return Spheroid(rev, refractive_index, axis_ratio, λ, ScattererInfo(T))
    elseif shape == SHAPE_CYLINDER
        return Cylinder(rev, refractive_index, axis_ratio, λ, ScattererInfo(T))
    elseif shape == SHAPE_CHEBYSHEV
        return Chebyshev(rev, refractive_index, axis_ratio, n, λ, ScattererInfo(T))
    end
end

function has_symmetric_plane(scatterer::AbstractScatterer)
    return typeof(scatterer) <: Spheroid ||
           typeof(scatterer) <: Cylinder ||
           (typeof(scatterer) <: Chebyshev && scatterer.n % 2 == 0)
end

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
    ngstart::Int64 = 0,
) where {T<:Real}
    kr = 2 * T(π) * scatterer.rev / scatterer.λ
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
    ngstart::Int64 = 0,
) where {T<:Real}
    kr = 2 * T(π) * scatterer.rev / scatterer.λ
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
            break
        else
            ngauss, nmax = nngauss, nnmax
        end
        GC.gc() # Trigger GC manually to avoid memory leak
    end

    @debug "Calculate T-Matrix for m = 0"
    T0, _ = tmatr0!(scatterer, ngauss, nmax)
    TT = [T0]
    for m in 1:nmax
        @debug "Calculate T-Matrix for m = $m"
        Tm, _ = tmatr!(scatterer, m, ngauss, nmax)
        push!(TT, Tm)
        GC.gc() # Trigger GC manually to avoid memory leak
    end

    @debug "Cross section" cross_section(TT, scatterer.λ)

    if T <: Arb
        @debug save_accuracy_info()
    end

    return TT
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
                    CN1 = CAL[n + dm, nn + dm] * fc
                    CN2 = CAL[n + dm, nn + dm] * fs
                    D11 = m^2 * dv1[n + dm] * dv01[nn + dm]
                    D12 = m * dv1[n + dm] * dv02[nn + dm]
                    D21 = m * dv2[n + dm] * dv01[nn + dm]
                    D22 = dv2[n + dm] * dv02[nn + dm]
                    VV += (T11 * D11 + T21 * D21 + T12 * D12 + T22 * D22) * CN1
                    VH += (T11 * D12 + T21 * D22 + T12 * D11 + T22 * D21) * CN2
                    HV -= (T11 * D21 + T21 * D11 + T12 * D22 + T22 * D12) * CN2
                    HH += (T11 * D22 + T21 * D12 + T12 * D21 + T22 * D11) * CN1
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
    hovenr(α₁, α₂, α₃, α₄, β₁, β₂)

    return α₁, α₂, α₃, α₄, β₁, β₂
end

@doc raw"""
```
calc_scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, Θ)
```

Calculate the scatterering matrix from the given expansion coefficients.

Parameters:

- `α₁`, `α₂`, `α₃`, `α₄`, `β₁`, `β₂`: The precalculated expansion coefficients.
- `Θ`: The scattering angle in degrees.
"""
function calc_scattering_matrix(
    α₁::AbstractVector{T},
    α₂::AbstractVector{T},
    α₃::AbstractVector{T},
    α₄::AbstractVector{T},
    β₁::AbstractVector{T},
    β₂::AbstractVector{T},
    Θ::Real,
) where {T<:Real}
    lmax = length(α₁) - 1
    Θ = Float64(Θ) / 180 * π

    F₁₁ = sum(α₁[l] * WignerD.wignerdjmn(l, 0, 0, Θ) for l in 0:lmax)
    F₂₂₊₃₃ = sum((α₂[l] + α₃[l]) * WignerD.wignerdjmn(l, 2, 2, Θ) for l in 2:lmax)
    F₂₂₋₃₃ = sum((α₂[l] - α₃[l]) * WignerD.wignerdjmn(l, 2, -2, Θ) for l in 2:lmax)
    F₂₂ = (F₂₂₊₃₃ + F₂₂₋₃₃) / 2
    F₃₃ = F₂₂₊₃₃ - F₂₂
    F₄₄ = sum(α₄[l] * WignerD.wignerdjmn(l, 0, 0, Θ) for l in 0:lmax)
    F₁₂ = -sum(β₁[l] * WignerD.wignerdjmn(l, 0, 2, Θ) for l in 2:lmax)
    F₃₄ = -sum(β₂[l] * WignerD.wignerdjmn(l, 0, 2, Θ) for l in 2:lmax)

    return F₁₁, F₂₂, F₃₃, F₄₄, F₁₂, F₃₄
end

function theta_split!(
    scatterer::AbstractScatterer{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
) where {T<:Real}
    if typeof(scatterer) <: Cylinder
        ng = ngauss ÷ 2
        ng1 = ng ÷ 2
        ng2 = ng - ng1
        x1, w1 = gausslegendre(T, ng1)
        x2, w2 = gausslegendre(T, ng2)
        xx = -cos(atan(scatterer.d_to_h))
        x[1:ng1] .= 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
        w[1:ng1] .= 0.5(xx + 1.0) .* w1
        x[(ng1 + 1):ng] .= -0.5xx .* x2 .+ 0.5xx
        w[(ng1 + 1):ng] .= -0.5xx .* w2
        x[(ng + 1):ngauss] .= (-1.0) .* x[ng:-1:1]
        w[(ng + 1):ngauss] .= w[ng:-1:1]
    else
        x0, w0 = gausslegendre(T, ngauss)
        x .= x0
        w .= w0
    end
end

function theta_split(scatterer::AbstractScatterer{T}, ngauss::Int64) where {T<:Real}
    x = zeros(ngauss)
    w = zeros(ngauss)
    theta_split!(scatterer, ngauss, x, w)
    return x, w
end

@doc raw"""
```
calc_r!(scatterer::AbstractScatterer{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer, in place.
"""
function calc_r!(
    scatterer::AbstractScatterer{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev

    if typeof(scatterer) <: Cylinder
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for cylinders")
        end

        e = scatterer.d_to_h
        h = rev * ∛(2 / (3e^2))
        d = h * e

        @simd for i in 1:(ngauss ÷ 2)
            cosθ = abs(x[i])
            sinθ = √(1 - cosθ^2)
            if h / cosθ < d / sinθ
                r[i] = h / cosθ
                dr[i] = h * sinθ / cosθ^2
            else
                r[i] = d / sinθ
                dr[i] = -d * cosθ / sinθ^2
            end
            r[ngauss + 1 - i] = r[i]
            dr[ngauss + 1 - i] = dr[i]
            dr[i] = -dr[i]
        end

    elseif typeof(scatterer) <: Spheroid
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for spheroids")
        end

        e = scatterer.a_to_c
        a = rev * ∛e

        @simd for i in 1:(ngauss ÷ 2)
            cosθ = x[i]
            sinθ = √(1.0 - cosθ^2)
            r[i] = a * √(1.0 / (e^2 * cosθ^2 + sinθ^2))
            r[ngauss + 1 - i] = r[i]
            dr[i] = r[i]^3 * cosθ * sinθ * (e^2 - 1.0) / a^2
            dr[ngauss + 1 - i] = -dr[i]
        end
    else
        @assert typeof(scatterer) <: Chebyshev
        e = scatterer.ε
        n = scatterer.n
        dn = float(n)

        a = 1.5e^2 * (4.0dn^2 - 2.0) / (4.0dn^2 - 1.0) + 1.0

        if n % 2 == 0
            a -= 3.0e * (1.0 + 0.25e^2) / (dn^2 - 1.0) + 0.25e^3 / (9.0dn^2 - 1.0)
        end

        r0 = rev / ∛a

        @simd for i in 1:ngauss
            xi = acos(x[i]) * n
            r[i] = r0 * (1.0 + e * cos(xi))
            dr[i] = -r0 * e * n * sin(xi)
        end
    end
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
            info.an = [T(n * (n + 1)) for n in 1:(info.ncap)]
            info.ann = [
                0.5 * √T((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1))) for n1 in 1:(info.ncap),
                n2 in 1:(info.ncap)
            ]
            info.sig = [i % 2 == 1 ? T(-1.0) : T(1.0) for i in 1:(info.ncap)]
        end

        if info.ngcap > ngcap
            info.x = zeros(T, info.ngcap)
            info.w = zeros(T, info.ngcap)
            info.s = zeros(T, info.ngcap)
            info.r = zeros(T, info.ngcap)
            info.dr = zeros(T, info.ngcap)
            info.kr1 = zeros(T, info.ngcap)
            info.kr_s1 = zeros(Complex{T}, info.ngcap)
        end

        if info.ncap > ncap
            info.j_tmp = zeros(T, 3info.ncap)
            info.j_s_tmp = zeros(Complex{T}, 3info.ncap)
            info.J11 = zeros(Complex{T}, info.ncap, info.ncap)
            info.J12 = zeros(Complex{T}, info.ncap, info.ncap)
            info.J21 = zeros(Complex{T}, info.ncap, info.ncap)
            info.J22 = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ11 = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ12 = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ21 = zeros(Complex{T}, info.ncap, info.ncap)
            info.RgJ22 = zeros(Complex{T}, info.ncap, info.ncap)
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
            info.jkr_s = zeros(Complex{T}, info.ngcap, info.ncap)
            info.djkr_s = zeros(Complex{T}, info.ngcap, info.ncap)
        end
    end

    # Need to recalculate `x`, `w`, `s`, `r`, `dr`, `kr1` and `kr_s1` only if `ngauss` changes.
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
        info.s[1:ngauss] .= 1 ./ (sin ∘ acos).(info.x[1:ngauss])
        info.kr1[1:ngauss] .= 1.0 ./ (k * info.r[1:ngauss])
        info.kr_s1[1:ngauss] .= 1.0 ./ (scatterer.m * k .* info.r[1:ngauss])
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.

    kr = k * info.r[1:ngauss]
    kr_s = scatterer.m * kr
    rmax = maximum(info.r[1:ngauss])
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(8.0 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb + 4.0 * ∛tb + 8.0 * √tb - nmax + 5))

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, nnmax1, view(info.jkr, i, :), view(info.djkr, i, :), info.j_tmp)
        sphericalbessely!(kr[i], nmax, view(info.ykr, i, :), view(info.dykr, i, :))
        sphericalbesselj!(kr_s[i], nmax, nnmax2, view(info.jkr_s, i, :), view(info.djkr_s, i, :), info.j_s_tmp)
    end

    view(info.hkr, 1:ngauss, 1:nmax) .= complex.(view(info.jkr, 1:ngauss, 1:nmax), view(info.ykr, 1:ngauss, 1:nmax))
    view(info.dhkr, 1:ngauss, 1:nmax) .= complex.(view(info.djkr, 1:ngauss, 1:nmax), view(info.dykr, 1:ngauss, 1:nmax))

    info.ngauss = ngauss
    info.nmax = nmax

    return
end

function update!(scatterer::AbstractScatterer{Arb}, ngauss::Int64, nmax::Int64)
    info = scatterer.info

    # No need to recalculate if both `ngauss` and `nmax` remains the same.
    if ngauss == info.ngauss && nmax == info.nmax
        return
    end

    # Need to recalculate `an`, `ann` and `sig` if `nmax` changes.
    if nmax != info.nmax
        info.an = ArbRefVector([Arb(n * (n + 1)) for n in 1:nmax])
        info.ann = ArbRefMatrix([
            √(Arb(2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1))) / 2 for n1 in 1:nmax, n2 in 1:nmax
        ])
        info.sig = ArbRefVector([i % 2 == 1 ? -1 : 1 for i in 1:nmax])
    end

    # Need to recalculate `x`, `w`, `s`, `r`, `dr`, `kr1` and `kr_s1` if `ngauss` changes.
    k = 2 * Arb(π) / scatterer.λ

    if ngauss != info.ngauss
        info.x = ArbRefVector(ngauss)
        info.w = ArbRefVector(ngauss)
        info.s = ArbRefVector(ngauss)
        info.r = ArbRefVector(ngauss)
        info.dr = ArbRefVector(ngauss)
        info.kr1 = ArbRefVector(ngauss)
        info.kr_s1 = AcbRefVector(ngauss)

        calc_r!(scatterer, ngauss, info.x, info.w, info.r, info.dr)
        info.s .= 1 ./ (sin ∘ acos).(info.x)
        info.kr1 .= 1 ./ (k * info.r)
        info.kr_s1 .= 1 ./ (scatterer.m * k .* info.r)
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.

    kr = ArbRefVector(ngauss)
    Arblib.mul!(kr, info.r, k)
    kr_s = AcbRefVector(ngauss)
    Arblib.mul!(kr_s, AcbRefVector(kr), scatterer.m)

    info.jkr = ArbRefMatrix(ngauss, nmax)
    info.djkr = ArbRefMatrix(ngauss, nmax)
    info.ykr = ArbRefMatrix(ngauss, nmax)
    info.dykr = ArbRefMatrix(ngauss, nmax)
    info.hkr = AcbRefMatrix(ngauss, nmax)
    info.dhkr = AcbRefMatrix(ngauss, nmax)
    info.jkr_s = AcbRefMatrix(ngauss, nmax)
    info.djkr_s = AcbRefMatrix(ngauss, nmax)

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, view(info.jkr, i, :), view(info.djkr, i, :))
        sphericalbessely!(kr[i], nmax, view(info.ykr, i, :), view(info.dykr, i, :))
        sphericalbesselj!(kr_s[i], nmax, view(info.jkr_s, i, :), view(info.djkr_s, i, :))
    end

    view(info.hkr, 1:ngauss, 1:nmax) .= complex.(view(info.jkr, 1:ngauss, 1:nmax), view(info.ykr, 1:ngauss, 1:nmax))
    view(info.dhkr, 1:ngauss, 1:nmax) .= complex.(view(info.djkr, 1:ngauss, 1:nmax), view(info.dykr, 1:ngauss, 1:nmax))

    info.ngauss = ngauss
    info.nmax = nmax
    info.ngcap = ngauss
    info.ncap = nmax

    return
end

function constant(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64) where {T<:Real}
    an = [Float64(n * (n + 1)) for n in 1:nmax]
    ann = [0.5 * √((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1))) for n1 in 1:nmax, n2 in 1:nmax]

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
    kr1 = 1.0 ./ kr
    kr_s = scatterer.m * kr
    kr_s1 = 1.0 ./ kr_s
    rmax = maximum(r)
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(8.0 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb + 4.0 * ∛tb + 8.0 * √tb - nmax + 5))

    jkr = zeros(T, ngauss, nmax)
    djkr = zeros(T, ngauss, nmax)
    ykr = zeros(T, ngauss, nmax)
    dykr = zeros(T, ngauss, nmax)
    jkr_s = zeros(Complex{T}, ngauss, nmax)
    djkr_s = zeros(Complex{T}, ngauss, nmax)
    j_tmp = zeros(T, nmax + nnmax1)
    j_s_tmp = zeros(Complex{T}, nmax + nnmax2)

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, nnmax1, view(jkr, i, :), view(djkr, i, :), j_tmp)
        sphericalbessely!(kr[i], nmax, view(ykr, i, :), view(dykr, i, :))
        sphericalbesselj!(kr_s[i], nmax, nnmax2, view(jkr_s, i, :), view(djkr_s, i, :), j_s_tmp)
    end

    return r, dr, kr1, kr_s1, jkr, djkr, ykr, dykr, jkr_s, djkr_s
end

function tmatr0!(scatterer::AbstractScatterer{T}, ngauss::Int64, nmax::Int64;) where {T<:Real}
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
        τ[ineg, :] .= view(τ, i, :) .* sig .* (-1.0)
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    w = view(info.w, 1:ngss)
    r = view(info.r, 1:ngss)
    dr = view(info.dr, 1:ngss)
    drr = dr ./ r
    kr1 = view(info.kr1, 1:ngss)
    kr_s1 = view(info.kr_s1, 1:ngss)
    wr2 = w .* r .* r

    jkr = view(info.jkr, 1:ngss, 1:nmax)
    djkr = view(info.djkr, 1:ngss, 1:nmax)
    hkr = view(info.hkr, 1:ngss, 1:nmax)
    dhkr = view(info.dhkr, 1:ngss, 1:nmax)
    jkr_s = view(info.jkr_s, 1:ngss, 1:nmax)
    djkr_s = view(info.djkr_s, 1:ngss, 1:nmax)

    J12 = view(info.J12, 1:nmax, 1:nmax)
    J21 = view(info.J21, 1:nmax, 1:nmax)
    RgJ12 = view(info.RgJ12, 1:nmax, 1:nmax)
    RgJ21 = view(info.RgJ21, 1:nmax, 1:nmax)
    fill!(J12, 0.0)
    fill!(J21, 0.0)
    fill!(RgJ12, 0.0)
    fill!(RgJ21, 0.0)

    Threads.@threads for n2 in 1:nmax
        for n1 in 1:nmax
            for i in 1:ngss
                τ₁τ₂ = τ[i, n1] * τ[i, n2]
                d₁τ₂ = d[i, n1] * τ[i, n2]
                d₂τ₁ = d[i, n2] * τ[i, n1]

                if !(sym && (n1 + n2) % 2 == 1)
                    J12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                    J21[n1, n2] +=
                        wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                    RgJ12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                    RgJ21[n1, n2] +=
                        wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
                end
            end
        end
    end

    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann

    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann

    k = 2 * T(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = k^2
    kk_s = k * k_s

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    Q = view(info.Q, 1:(2nmax), 1:(2nmax))
    Q11 = view(Q, 1:nmax, 1:nmax)
    Q22 = view(Q, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    RgQ = view(info.RgQ, 1:(2nmax), 1:(2nmax))
    RgQ11 = view(RgQ, 1:nmax, 1:nmax)
    RgQ22 = view(RgQ, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    fill!(Q, zero(eltype(Q)))
    fill!(RgQ, zero(eltype(RgQ)))

    Q11 .= kk_s * J21 .+ kk * J12
    Q22 .= kk_s * J12 .+ kk * J21

    RgQ11 .= kk_s * RgJ21 .+ kk * RgJ12
    RgQ22 .= kk_s * RgJ12 .+ kk * RgJ21

    T0 = -RgQ * inv(Q)

    return T0, Q, RgQ
end

function tmatr!(scatterer::AbstractScatterer{T}, m::Int64, ngauss::Int64, nmax::Int64;) where {T<:Real}
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
    w = view(info.w, 1:ngss)
    r = view(info.r, 1:ngss)
    dr = view(info.dr, 1:ngss)
    drr = dr ./ r
    kr1 = view(info.kr1, 1:ngss)
    kr_s1 = view(info.kr_s1, 1:ngss)
    wr2 = w .* r .* r

    jkr = view(info.jkr, 1:ngss, 1:nmax)
    djkr = view(info.djkr, 1:ngss, 1:nmax)
    hkr = view(info.hkr, 1:ngss, 1:nmax)
    dhkr = view(info.dhkr, 1:ngss, 1:nmax)
    jkr_s = view(info.jkr_s, 1:ngss, 1:nmax)
    djkr_s = view(info.djkr_s, 1:ngss, 1:nmax)

    J11 = view(info.J11, mm:nmax, mm:nmax)
    J12 = view(info.J12, mm:nmax, mm:nmax)
    J21 = view(info.J21, mm:nmax, mm:nmax)
    J22 = view(info.J22, mm:nmax, mm:nmax)
    RgJ11 = view(info.RgJ11, mm:nmax, mm:nmax)
    RgJ12 = view(info.RgJ12, mm:nmax, mm:nmax)
    RgJ21 = view(info.RgJ21, mm:nmax, mm:nmax)
    RgJ22 = view(info.RgJ22, mm:nmax, mm:nmax)
    fill!(J11, zero(eltype(J11)))
    fill!(J12, zero(eltype(J12)))
    fill!(J21, zero(eltype(J21)))
    fill!(J22, zero(eltype(J22)))
    fill!(RgJ11, zero(eltype(RgJ11)))
    fill!(RgJ12, zero(eltype(RgJ12)))
    fill!(RgJ21, zero(eltype(RgJ21)))
    fill!(RgJ22, zero(eltype(RgJ22)))

    OffsetJ11 = OffsetArray(J11, mm:nmax, mm:nmax)
    OffsetJ12 = OffsetArray(J12, mm:nmax, mm:nmax)
    OffsetJ21 = OffsetArray(J21, mm:nmax, mm:nmax)
    OffsetJ22 = OffsetArray(J22, mm:nmax, mm:nmax)
    OffsetRgJ11 = OffsetArray(RgJ11, mm:nmax, mm:nmax)
    OffsetRgJ12 = OffsetArray(RgJ12, mm:nmax, mm:nmax)
    OffsetRgJ21 = OffsetArray(RgJ21, mm:nmax, mm:nmax)
    OffsetRgJ22 = OffsetArray(RgJ22, mm:nmax, mm:nmax)

    Threads.@threads for n2 in mm:nmax
        for n1 in mm:nmax
            for i in 1:ngss
                if !(sym && (n1 + n2) % 2 == 0)
                    pττp = p[i, n1] * τ[i, n2] + p[i, n2] * τ[i, n1]
                    p₁d₂ = p[i, n1] * d[i, n2]

                    OffsetJ11[n1, n2] += wr2[i] * hkr[i, n1] * jkr_s[i, n2] * pττp

                    OffsetJ22[n1, n2] +=
                        wr2[i] * (
                            dhkr[i, n1] * djkr_s[i, n2] * pττp +
                            drr[i] *
                            (
                                an[n1] * hkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * dhkr[i, n1]
                            ) *
                            p₁d₂
                        )

                    OffsetRgJ11[n1, n2] += wr2[i] * jkr[i, n1] * jkr_s[i, n2] * pττp

                    OffsetRgJ22[n1, n2] +=
                        wr2[i] * (
                            djkr[i, n1] * djkr_s[i, n2] * pττp +
                            drr[i] *
                            (
                                an[n1] * jkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * djkr[i, n1]
                            ) *
                            p₁d₂
                        )
                end

                if !(sym && (n1 + n2) % 2 == 1)
                    ppττ = p[i, n1] * p[i, n2] + τ[i, n1] * τ[i, n2]
                    d₁τ₂ = d[i, n1] * τ[i, n2]
                    d₂τ₁ = d[i, n2] * τ[i, n1]

                    OffsetJ12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * ppττ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                    OffsetJ21[n1, n2] +=
                        wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                    OffsetRgJ12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * ppττ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                    OffsetRgJ21[n1, n2] +=
                        wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
                end
            end
        end
    end

    ann = view(ann, mm:nmax, mm:nmax)
    J11 .*= -ann
    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann
    J22 .*= -ann

    RgJ11 .*= -ann
    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann
    RgJ22 .*= -ann

    k = 2 * T(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = k^2
    kk_s = k * k_s

    nm = nmax - mm + 1
    Q = view(info.Q, 1:(2nm), 1:(2nm))
    Q11 = view(Q, 1:nm, 1:nm)
    Q12 = view(Q, 1:nm, (nm + 1):(2nm))
    Q21 = view(Q, (nm + 1):(2nm), 1:nm)
    Q22 = view(Q, (nm + 1):(2nm), (nm + 1):(2nm))
    RgQ = view(info.RgQ, 1:(2nm), 1:(2nm))
    RgQ11 = view(RgQ, 1:nm, 1:nm)
    RgQ12 = view(RgQ, 1:nm, (nm + 1):(2nm))
    RgQ21 = view(RgQ, (nm + 1):(2nm), 1:nm)
    RgQ22 = view(RgQ, (nm + 1):(2nm), (nm + 1):(2nm))
    fill!(Q, zero(eltype(Q)))
    fill!(RgQ, zero(eltype(RgQ)))

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    Q11 .= kk_s * J21 .+ kk * J12
    Q12 .= kk_s * J11 .+ kk * J22
    Q21 .= kk_s * J22 .+ kk * J11
    Q22 .= kk_s * J12 .+ kk * J21

    RgQ11 .= kk_s * RgJ21 .+ kk * RgJ12
    RgQ12 .= kk_s * RgJ11 .+ kk * RgJ22
    RgQ21 .= kk_s * RgJ22 .+ kk * RgJ11
    RgQ22 .= kk_s * RgJ12 .+ kk * RgJ21

    Tm = -RgQ * inv(Q)

    return Tm, Q, RgQ
end

function tmatr0!(scatterer::AbstractScatterer{Arb}, ngauss::Int64, nmax::Int64)
    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    an = info.an
    ann = info.ann
    sig = info.sig
    x = info.x

    d = [ArbRefVector(nmax) for _ in 1:ngauss]
    τ = [ArbRefVector(nmax) for _ in 1:ngauss]
    for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, 0, x[i], d[i], τ[i])
        @. d[ineg] = d[i] * sig
        @. τ[ineg] = τ[i] * sig * (-1)
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    w = info.w
    r = info.r
    dr = info.dr
    drr = ArbRefVector(dr ./ r)
    kr1 = info.kr1
    kr_s1 = info.kr_s1
    wr2 = ArbRefVector(w .* r .* r)

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkr_s = info.jkr_s
    djkr_s = info.djkr_s

    J12 = AcbRefMatrix(nmax, nmax)
    J21 = AcbRefMatrix(nmax, nmax)
    RgJ12 = AcbRefMatrix(nmax, nmax)
    RgJ21 = AcbRefMatrix(nmax, nmax)

    Threads.@threads for n2 in 1:nmax
        for n1 in 1:nmax
            for i in 1:ngss
                τ₁τ₂ = τ[i][n1] * τ[i][n2]
                d₁τ₂ = d[i][n1] * τ[i][n2]
                d₂τ₁ = d[i][n2] * τ[i][n1]

                if !(sym && (n1 + n2) % 2 == 1)
                    J12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                    J21[n1, n2] +=
                        wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                    RgJ12[n1, n2] +=
                        wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                    RgJ21[n1, n2] +=
                        wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
                end
            end
        end
    end

    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann

    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann

    k = 2 * Arb(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = Acb(k^2)
    kk_s = k * k_s

    Q11 = AcbRefMatrix(nmax, nmax)
    Q12 = AcbRefMatrix(nmax, nmax)
    Q21 = AcbRefMatrix(nmax, nmax)
    Q22 = AcbRefMatrix(nmax, nmax)
    RgQ11 = AcbRefMatrix(nmax, nmax)
    RgQ12 = AcbRefMatrix(nmax, nmax)
    RgQ21 = AcbRefMatrix(nmax, nmax)
    RgQ22 = AcbRefMatrix(nmax, nmax)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nmax, nmax)
    Arblib.mul!(Q11, J21, kk_s)
    Arblib.mul!(tmp, J12, kk)
    Arblib.add!(Q11, Q11, tmp)

    Arblib.mul!(Q22, J12, kk_s)
    Arblib.mul!(tmp, J21, kk)
    Arblib.add!(Q22, Q22, tmp)

    Arblib.mul!(RgQ11, RgJ21, kk_s)
    Arblib.mul!(tmp, RgJ12, kk)
    Arblib.add!(RgQ11, RgQ11, tmp)

    Arblib.mul!(RgQ22, RgJ12, kk_s)
    Arblib.mul!(tmp, RgJ21, kk)
    Arblib.add!(RgQ22, RgQ22, tmp)

    Q = vcat(hcat(Q11, Q12), hcat(Q21, Q22))
    RgQ = vcat(hcat(RgQ11, RgQ12), hcat(RgQ21, RgQ22))

    IQ = similar(Q)
    if ARB_APPROX_INV[]
        ret = Arblib.approx_inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        T0 = -RgQ * IQ
    else
        ret = Arblib.inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        T0 = -RgQ * IQ

        if COLLECT_ACCURACY_INFO[]
            @debug begin
                push!(ACCURACY_INFO[], ("d", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(d)))
                push!(ACCURACY_INFO[], ("τ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(τ)))
                push!(ACCURACY_INFO[], ("J12", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J12)))
                push!(ACCURACY_INFO[], ("J21", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J21)))
                push!(ACCURACY_INFO[], ("RgJ12", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ12)))
                push!(ACCURACY_INFO[], ("RgJ21", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ21)))
                push!(ACCURACY_INFO[], ("RgQ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgQ)))
                push!(ACCURACY_INFO[], ("Q", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(Q)))
                push!(ACCURACY_INFO[], ("inv(Q)", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(IQ)))
                push!(ACCURACY_INFO[], ("T", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(T0)))
                "Collecting accuracy info..."
            end
        end

        @debug "Accuracy of T is $(rel_accuracy_bits(T0))"
    end

    return T0, Q, RgQ
end

function tmatr!(scatterer::AbstractScatterer{Arb}, m::Int64, ngauss::Int64, nmax::Int64;)
    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    mm = max(m, 1)
    an = info.an
    ann = info.ann
    sig = info.sig
    x = info.x
    s = info.s

    d = [ArbRefVector(nmax) for _ in 1:ngauss]
    p = [ArbRefVector(nmax) for _ in 1:ngauss]
    τ = [ArbRefVector(nmax) for _ in 1:ngauss]

    for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, m, x[i], d[i], τ[i])
        @. p[i] = d[i] * s[i] * m
        @. d[ineg] = d[i] * sig
        @. τ[ineg] = τ[i] * sig * (-1)
        @. p[ineg] = p[i] * sig
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    w = info.w
    r = info.r
    dr = info.dr
    drr = ArbRefVector(dr ./ r)
    kr1 = info.kr1
    kr_s1 = info.kr_s1
    wr2 = ArbRefVector(w .* r .* r)

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkr_s = info.jkr_s
    djkr_s = info.djkr_s

    nm = nmax - mm + 1
    J11 = AcbRefMatrix(nm, nm)
    J12 = AcbRefMatrix(nm, nm)
    J21 = AcbRefMatrix(nm, nm)
    J22 = AcbRefMatrix(nm, nm)
    RgJ11 = AcbRefMatrix(nm, nm)
    RgJ12 = AcbRefMatrix(nm, nm)
    RgJ21 = AcbRefMatrix(nm, nm)
    RgJ22 = AcbRefMatrix(nm, nm)

    Threads.@threads for n2 in mm:nmax
        nn2 = n2 - mm + 1
        for n1 in mm:nmax
            nn1 = n1 - mm + 1
            for i in 1:ngss
                if !(sym && (n1 + n2) % 2 == 0)
                    pττp = p[i][n1] * τ[i][n2] + p[i][n2] * τ[i][n1]
                    p₁d₂ = p[i][n1] * d[i][n2]

                    J11[nn1, nn2] += wr2[i] * hkr[i, n1] * jkr_s[i, n2] * pττp

                    J22[nn1, nn2] +=
                        wr2[i] * (
                            dhkr[i, n1] * djkr_s[i, n2] * pττp +
                            drr[i] *
                            (
                                an[n1] * hkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * dhkr[i, n1]
                            ) *
                            p₁d₂
                        )

                    RgJ11[nn1, nn2] += wr2[i] * jkr[i, n1] * jkr_s[i, n2] * pττp

                    RgJ22[nn1, nn2] +=
                        wr2[i] * (
                            djkr[i, n1] * djkr_s[i, n2] * pττp +
                            drr[i] *
                            (
                                an[n1] * jkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * djkr[i, n1]
                            ) *
                            p₁d₂
                        )
                end

                if !(sym && (n1 + n2) % 2 == 1)
                    ppττ = p[i][n1] * p[i][n2] + τ[i][n1] * τ[i][n2]
                    d₁τ₂ = d[i][n1] * τ[i][n2]
                    d₂τ₁ = d[i][n2] * τ[i][n1]

                    J12[nn1, nn2] +=
                        wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * ppττ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                    J21[nn1, nn2] +=
                        wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                    RgJ12[nn1, nn2] +=
                        wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * ppττ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                    RgJ21[nn1, nn2] +=
                        wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
                end
            end
        end
    end

    ann = view(ann, mm:nmax, mm:nmax)
    J11 .*= -ann
    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann
    J22 .*= -ann

    RgJ11 .*= -ann
    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann
    RgJ22 .*= -ann

    k = 2 * Arb(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = Acb(k^2)
    kk_s = k * k_s

    Q11 = AcbRefMatrix(nm, nm)
    Q12 = AcbRefMatrix(nm, nm)
    Q21 = AcbRefMatrix(nm, nm)
    Q22 = AcbRefMatrix(nm, nm)
    RgQ11 = AcbRefMatrix(nm, nm)
    RgQ12 = AcbRefMatrix(nm, nm)
    RgQ21 = AcbRefMatrix(nm, nm)
    RgQ22 = AcbRefMatrix(nm, nm)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nm, nm)
    Arblib.mul!(Q11, J21, kk_s)
    Arblib.mul!(tmp, J12, kk)
    Arblib.add!(Q11, Q11, tmp)

    Arblib.mul!(Q12, J11, kk_s)
    Arblib.mul!(tmp, J22, kk)
    Arblib.add!(Q12, Q12, tmp)

    Arblib.mul!(Q21, J22, kk_s)
    Arblib.mul!(tmp, J11, kk)
    Arblib.add!(Q21, Q21, tmp)

    Arblib.mul!(Q22, J12, kk_s)
    Arblib.mul!(tmp, J21, kk)
    Arblib.add!(Q22, Q22, tmp)

    Arblib.mul!(RgQ11, RgJ21, kk_s)
    Arblib.mul!(tmp, RgJ12, kk)
    Arblib.add!(RgQ11, RgQ11, tmp)

    Arblib.mul!(RgQ12, RgJ11, kk_s)
    Arblib.mul!(tmp, RgJ22, kk)
    Arblib.add!(RgQ12, RgQ12, tmp)

    Arblib.mul!(RgQ21, RgJ22, kk_s)
    Arblib.mul!(tmp, RgJ11, kk)
    Arblib.add!(RgQ21, RgQ21, tmp)

    Arblib.mul!(RgQ22, RgJ12, kk_s)
    Arblib.mul!(tmp, RgJ21, kk)
    Arblib.add!(RgQ22, RgQ22, tmp)

    Q = vcat(hcat(Q11, Q12), hcat(Q21, Q22))
    RgQ = vcat(hcat(RgQ11, RgQ12), hcat(RgQ21, RgQ22))

    IQ = similar(Q)
    if ARB_APPROX_INV[]
        ret = Arblib.approx_inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        Tm = -RgQ * IQ
    else
        ret = Arblib.inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        Tm = -RgQ * IQ

        if COLLECT_ACCURACY_INFO[]
            @debug begin
                push!(ACCURACY_INFO[], ("d", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(d)))
                push!(ACCURACY_INFO[], ("p", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(p)))
                push!(ACCURACY_INFO[], ("τ", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(τ)))
                push!(ACCURACY_INFO[], ("J11", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J11)))
                push!(ACCURACY_INFO[], ("J12", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J12)))
                push!(ACCURACY_INFO[], ("J21", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J21)))
                push!(ACCURACY_INFO[], ("J22", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J22)))
                push!(ACCURACY_INFO[], ("RgJ11", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ11)))
                push!(ACCURACY_INFO[], ("RgJ12", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ12)))
                push!(ACCURACY_INFO[], ("RgJ21", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ21)))
                push!(ACCURACY_INFO[], ("RgJ22", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ22)))
                push!(ACCURACY_INFO[], ("RgQ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgQ)))
                push!(ACCURACY_INFO[], ("Q", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(Q)))
                push!(ACCURACY_INFO[], ("inv(Q)", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(IQ)))
                push!(ACCURACY_INFO[], ("T", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(Tm)))
                "Collecting accuracy info..."
            end
        end

        @debug "Accuracy of T is $(rel_accuracy_bits(Tm))"
    end

    return Tm, Q, RgQ
end
