@enum Shape SHAPE_SPHEROID SHAPE_CYLINDER SHAPE_CHEBYSHEV
@enum RadiusType RADIUS_EQUAL_VOLUME RADIUS_EQUAL_AREA RADIUS_MAXIMUM

const NPN1 = 100
const NPNG1 = 500

mutable struct ScattererInfo{T<:Real}
    nmax::Int64
    ngauss::Int64
    ncap::Int64
    ngcap::Int64
    an::Array{T,1} # (nmax,)
    ann::Array{T,2} # (nmax, nmax)
    sig::Array{T,1} # (nmax,)
    x::Array{T,1} # (ngauss,)
    w::Array{T,1} # (ngauss,)
    s::Array{T,1} # (ngauss,)
    r::Array{T,1} # (ngauss,)
    dr::Array{T,1} # (ngauss,)
    kr1::Array{T,1} # (ngauss,)
    kr_s1::Array{Complex{T},1} # (ngauss,)
    d::Array{T,2} # (ngauss, nmax)
    τ::Array{T,2} # (ngauss, nmax)
    p::Array{T,2} # (ngauss, nmax)
    jkr::Array{T,2} # (ngauss, nmax)
    djkr::Array{T,2} # (ngauss, nmax)
    j_tmp::Array{T,1} # (2ncap,)
    ykr::Array{T,2} # (ngauss, nmax)
    dykr::Array{T,2} # (ngauss, nmax)
    hkr::Array{Complex{T},2} # (ngauss, nmax)
    dhkr::Array{Complex{T},2} # (ngauss, nmax)
    jkr_s::Array{Complex{T},2} # (ngauss, nmax)
    djkr_s::Array{Complex{T},2} # (ngauss, nmax)
    j_s_tmp::Array{Complex{T},1} # (2ncap,)
    J11::Array{Complex{T},2} # (nmax, nmax)
    J12::Array{Complex{T},2} # (nmax, nmax)
    J21::Array{Complex{T},2} # (nmax, nmax)
    J22::Array{Complex{T},2} # (nmax, nmax)
    RgJ11::Array{Complex{T},2} # (nmax, nmax)
    RgJ12::Array{Complex{T},2} # (nmax, nmax)
    RgJ21::Array{Complex{T},2} # (nmax, nmax)
    RgJ22::Array{Complex{T},2} # (nmax, nmax)
    Q::Array{Complex{T},2} # (2nmax, 2nmax)
    RgQ::Array{Complex{T},2} # (2nmax, 2nmax)
end

function ScattererInfo(T)
    return ScattererInfo(
        0,
        0,
        NPN1,
        NPNG1,
        [T(n * (n + 1)) for n in 1:NPN1],
        [T(0.5 * √((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1)))) for n1 in 1:NPN1, n2 in 1:NPN1],
        [i % 2 == 1 ? T(-1.0) : T(1.0) for i in 1:NPN1],
        zeros(T, NPNG1),
        zeros(T, NPNG1),
        zeros(T, NPNG1),
        zeros(T, NPNG1),
        zeros(T, NPNG1),
        zeros(T, NPNG1),
        zeros(Complex{T}, NPNG1),
        zeros(T, NPNG1, NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(T, 2NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(T, NPNG1, NPN1),
        zeros(Complex{T}, NPNG1, NPN1),
        zeros(Complex{T}, NPNG1, NPN1),
        zeros(Complex{T}, NPNG1, NPN1),
        zeros(Complex{T}, NPNG1, NPN1),
        zeros(Complex{T}, 2NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, NPN1, NPN1),
        zeros(Complex{T}, 2NPN1, 2NPN1),
        zeros(Complex{T}, 2NPN1, 2NPN1),
    )
end

abstract type AbstractScatterer end

const CHEBYSHEV_DEFAULT_GAUSSIAN_POINTS = 60

@doc raw"""
A spheroid scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `a_to_c`: The ratio $a/c$ of the horizontal to rotational axes.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Spheroid{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    a_to_c::T
    λ::T
    info::ScattererInfo{T}
end

@doc raw"""
A cylinder scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `d_to_h`: The diameter-to-height ratio $D/H$.
- `λ`: The wavelength of the incident wave.
"""
struct Cylinder{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    d_to_h::T
    λ::T
    info::ScattererInfo{T}
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
"""
struct Chebyshev{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    ε::T
    n::Int64
    λ::T
    info::ScattererInfo{T}
end

@doc raw"""
Scatterer constructor with named parameters.

Parameters:

- `r`: The equivalent radius.
- `shape`: The particle shape. Possible values are `SHAPE_SPHEROID` (default), `SHAPE_CYLINDER` and `SHAPE_CHEBYSHEV`.
- `axis_ratio`: For spheroids, it is the ratio $a/b$ of the horizontal to rotational axes. For cylinders, it is the diameter-to-length ratio $D/L$. For Chebyshev particles, it is the deformation parameter $\varepsilon$.
- `radius_type`: The type of the equivalent radius. Possible values are `RADIUS_EQUAL_VOLUME` (default), `RADIUS_EQUAL_AREA` and `RADIUS_MAXIMUM`. All radius types will be transformed into the equivalent volume radius.
- `refractive_index`: The complex refractive index.
- `ngauss`: Number of points for Gaussian integration. Required only for Chebyshev particles.
- `n`: Degree of the Chebyshev polynomial. Required only for Chebyshev particles.

"""
function Scatterer(;
    r::T,
    shape::Shape = SHAPE_SPHEROID,
    axis_ratio::T,
    radius_type::RadiusType = RADIUS_EQUAL_VOLUME,
    refractive_index::Complex{T} = 1.0,
    n::Int64 = 2,
    ngauss::Int64 = CHEBYSHEV_DEFAULT_GAUSSIAN_POINTS,
    λ::T = 1.0,
) where {T<:Real}
    if shape == SHAPE_CHEBYSHEV && (abs(axis_ratio) < 0.0 || abs(axis_ratio) >= 1.0)
        error("Constraint violated: Chebyshev particles should have 0≤|ε|<1.")
    end

    if radius_type == RADIUS_EQUAL_VOLUME
        rev = r
    elseif radius_type == RADIUS_EQUAL_AREA
        if shape == SHAPE_SPHEROID
            d = axis_ratio
            if abs(d - 1.0) < eps(d)
                ratio = 1.0
            elseif d > 1.0
                e = √(1.0 - 1.0 / d^2)
                ratio = 1.0 / √(0.25(2.0d^(2.0 / 3.0) + d^(-4.0 / 3.0) * log((1.0 + e) / (1.0 - e)) / e))
            elseif d < 1.0
                e = √(1.0 - d^2)
                ratio = 1.0 / √(0.5(d^(2.0 / 3.0) + d^(-1.0 / 3.0) * asin(e) / e))
            end
        elseif shape == SHAPE_CYLINDER
            e = axis_ratio
            ratio = ∛(1.5 / e) / √((e + 2.0) / 2.0e)
        elseif shape == SHAPE_CHEBYSHEV
            e = axis_ratio
            en = e * n
            x, w = gausslegendre(ngauss)
            s = 0.0
            v = 0.0
            @simd for i in 1:ngauss
                θ = acos(x[i])
                nθ = n * θ
                sinθ = sin(θ)
                sinnθ = sin(nθ)
                cosnθ = cos(nθ)
                a = 1.0 + e * cosnθ
                ens = en * sinnθ
                s += w[i] * a * √(a^2 + ens^2)
                v += w[i] * (sinθ * a + x[i] * ens) * sinθ * a^2
            end
            rs = √(0.5s)
            rv = ∛(0.75v)
            ratio = rv / rs
        end

        rev = r * ratio
    elseif radius_type == RADIUS_MAXIMUM
        if shape == SHAPE_SPHEROID
            rev = axis_ratio > 1.0 ? r / ∛axis_ratio : r * axis_ratio^(2.0 / 3.0)
        elseif shape == SHAPE_CYLINDER
            rev = axis_ratio > 1.0 ? r * ∛(1.5 / axis_ratio) : r * ∛(1.5axis_ratio^2)
        elseif shape == SHAPE_CHEBYSHEV
            e = axis_ratio
            en = e * n
            x, w = gausslegendre(ngauss)
            θ = acos.(x)
            nθ = n * θ
            sinθ = sin.(θ)
            sinnθ = sin.(nθ)
            cosnθ = cos.(nθ)
            a = e * cosnθ .+ 1.0
            ens = en * sinnθ
            v = sum(w .* (sinθ .* a + x .* ens) .* sinθ * a .^ 2)
            rv = ∛(0.75v)
            rev = r / (1.0 + e) * rv
        end
    end

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
calc_tmatrix(scatterer::Scatterer, accuracy::Float64=0.001)
```

Calculate the T-Matrix of the scatterer.

Parameters:

- `scatterer`: The scatterer.
- `accuracy`: The intended threshold for convergence. Default is 0.001.

"""
function calc_tmatrix(scatterer::AbstractScatterer, accuracy::Float64 = 0.0001)
    kr = 2π * scatterer.rev / scatterer.λ
    nstart = max(4, Int64(floor(kr + 4.05 * ∛kr)))
    ngstart = nstart * 4
    nmax = nstart
    ngauss = ngstart
    Δ = 1.0

    Qext0 = 0.0
    Qsca0 = 0.0
    while Δ > accuracy
        Qext = 0.0
        Qsca = 0.0

        T0, _ = tmatr0!(scatterer, ngauss, nmax)
        for n in 1:nmax
            Qsca += (2n + 1) * real(T0[n, n] * T0[n, n]' + T0[n + nmax, n + nmax] * T0[n + nmax, n + nmax]')
            Qext += (2n + 1) * (real(T0[n, n] + real(T0[n + nmax, n + nmax])))
        end

        ΔQsca = abs((Qsca0 - Qsca) / Qsca)
        ΔQext = abs((Qext0 - Qext) / Qext)
        Δ = max(ΔQsca, ΔQext)

        @debug "nmax iteration" nmax ΔQsca ΔQext

        if Δ <= accuracy
            break
        end

        Qsca0 = Qsca
        Qext0 = Qext
        nmax += 1
        ngauss += 4
    end

    Δ = 1.0
    Qext0 = 0.0
    Qsca0 = 0.0
    while Δ > accuracy
        Qext = 0.0
        Qsca = 0.0

        T0, _ = tmatr0!(scatterer, ngauss, nmax)
        for n in 1:nmax
            Qsca += (2n + 1) * real(T0[n, n] * T0[n, n]' + T0[n + nmax, n + nmax] * T0[n + nmax, n + nmax]')
            Qext += (2n + 1) * (real(T0[n, n]) + real(T0[n + nmax, n + nmax]))
        end

        ΔQsca = abs((Qsca0 - Qsca) / Qsca)
        ΔQext = abs((Qext0 - Qext) / Qext)
        Δ = max(ΔQsca, ΔQext)

        @debug "ngauss iteration" ngauss ΔQsca ΔQext

        Qsca0 = Qsca
        Qext0 = Qext

        if Δ <= accuracy
            break
        end

        ngauss += 2
    end

    @debug "Convergence reached" nmax ngauss

    T0, _ = tmatr0!(scatterer, ngauss, nmax)
    T = [T0]
    for mm in 1:nmax
        Tmm, _ = tmatr!(scatterer, mm, ngauss, nmax)
        push!(T, Tmm)
    end

    @debug "Cross section" cross_section(T, scatterer.λ)

    return T
end

@doc raw"""
```
calc_amplitude(scatterer::Scatterer, tmatrix::Union{Array{Float64,2},Nothing}, ϑ_i::Float64, φ_i::Float64, ϑ_s::Float64, φ_s::Float64)
```

Calculate the amplitude matrix and the phase matrix, given the scatterer and the geometry of the incident and the scattered beam. Use pre-computed T-Matrix when possible.

Parameters:

- `scatterer`: The scatterer.
- `α, β`: The Euler angle.
- `ϑ_i`: The zenith angle of the incident beam.
- `ϑ_s`: The zenith angle of the scattered beam.
- `φ_i`: The azimuth angle of the indicent beam.
- `φ_s`: The azimuth angle of the scatterer beam.
- `tmatrix`: The pre-computed T-Matrix of the scatterer, or nothing if there is no pre-computation.

> All the angles here are input in degrees.
"""
function calc_amplitude(
    scatterer::AbstractScatterer,
    α::T,
    β::T,
    ϑ_i::T,
    ϑ_s::T,
    φ_i::T,
    φ_s::T,
    tmatrix::Union{Vector{Array{Complex{T},2}},Nothing} = nothing,
) where {T<:Real}
    # Validate the input angles
    @assert 0.0 <= α <= 360.0 &&
            0.0 <= β <= 180.0 &&
            0.0 <= ϑ_i <= 180.0 &&
            0.0 <= ϑ_s <= 180.0 &&
            0.0 <= φ_i <= 360.0 &&
            0.0 <= φ_s <= 360.0

    if tmatrix === nothing
        tmatrix = calc_tmatrix(scatterer)
    end

    nmax = length(tmatrix) - 1

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
                T11 = tmatrix[m + 1][n, nn]
                T22 = tmatrix[m + 1][n + nm, nn + nm]
                if m == 0
                    CN = CAL[n + dm, nn + dm] * dv2[n + dm] * dv02[nn + dm]
                    VV += CN * T22
                    HH += CN * T11
                else
                    T12 = tmatrix[m + 1][n, nn + nm]
                    T21 = tmatrix[m + 1][n + nm, nn]
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

    DK = 2π / scatterer.λ
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

Calculate the phase matrix using the given amplitude matrix $\mathbf{S}$.
"""
function calc_phase(S::Array{Complex{T},2}) where {T<:Real}
    @assert size(S) == (2, 2)

    Z = zeros(Complex{T}, 4, 4)
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

calc_Z = calc_phase

function calc_SZ(
    scatterer::AbstractScatterer,
    α::T,
    β::T,
    ϑ_i::T,
    ϑ_s::T,
    φ_i::T,
    φ_s::T,
    tmatrix::Union{Vector{Array{Complex{T},2}},Nothing} = nothing,
) where {T<:Real}
    S = calc_S(scatterer, α, β, ϑ_i, ϑ_s, φ_i, φ_s, tmatrix)
    Z = calc_Z(S)

    return S, Z
end

@doc raw"""
```
calc_r!(scatterer::Scatterer, ngauss::Int64, r::AbstractArray{Float64}, dr::AbstractArray{Float64})
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer.
"""
function calc_r!(scatterer::AbstractScatterer, ngauss::Int64, r::AbstractArray{Float64}, dr::AbstractArray{Float64})
    rev = scatterer.rev

    if typeof(scatterer) <: Cylinder
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for cylinders")
        end

        # For cylinders, r(theta) is not continuous, thus must be handled separately.
        ng1 = ngauss ÷ 4
        ng2 = ngauss ÷ 2 - ng1
        x1, _ = gausslegendre(ng1)
        x2, _ = gausslegendre(ng2)
        x = zeros(ngauss)
        xx = -cos(atan(scatterer.d_to_h))
        x[1:ng1] = 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
        x[(ng1 + 1):(ng1 + ng2)] = -0.5xx .* x2 .+ 0.5xx
        x[(ngauss ÷ 2 + 1):ngauss] = -x[(ngauss ÷ 2):-1:1]

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

        x, _ = gausslegendre(ngauss)

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

        x, _ = gausslegendre(ngauss)

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
function calc_r(scatterer::AbstractScatterer, ngauss::Int64)
    rev = scatterer.rev
    r = zeros(ngauss)
    dr = zeros(ngauss)
    calc_r!(scatterer, ngauss, r, dr)
    return r, dr
end

function update!(scatterer::AbstractScatterer, ngauss::Int64, nmax::Int64)
    info = scatterer.info

    # No need to recalculate if both `ngauss` and `nmax` remains the same.
    if ngauss == info.ngauss && nmax == info.nmax
        return
    end

    T = typeof(scatterer.rev)

    # Need to enlarge the arrays if either `ngauss` or `nmax` exceeds the current capacity.
    if ngauss > info.ngcap || nmax > info.ncap
        ncap = info.ncap
        ngcap = info.ngcap
        info.ncap = max(ncap, nmax * 2)
        info.ngcap = max(ngcap, ngauss * 2)

        if info.ncap > ncap
            info.an = [T(n * (n + 1)) for n in 1:(info.ncap)]
            info.ann = [
                T(0.5 * √((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1)))) for n1 in 1:(info.ncap),
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
            info.j_tmp = zeros(T, 2info.ncap)
            info.j_s_tmp = zeros(Complex{T}, 2info.ncap)
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
    k = 2π / scatterer.λ

    if ngauss != info.ngauss
        if typeof(scatterer) <: Cylinder
            ng = ngauss ÷ 2
            ng1 = ng ÷ 2
            ng2 = ng - ng1
            x1, w1 = gausslegendre(ng1)
            x2, w2 = gausslegendre(ng2)
            xx = -cos(atan(scatterer.d_to_h))
            info.x[1:ng1] = 0.5(xx + 1.0) * x1 .+ 0.5(xx - 1.0)
            info.w[1:ng1] = 0.5(xx + 1.0) * w1
            info.x[(ng1 + 1):ng] = -0.5xx * x2 .+ 0.5xx
            info.w[(ng1 + 1):ng] = -0.5xx * w2
            info.x[(ng + 1):ngauss] = -info.x[ng:-1:1]
            info.w[(ng + 1):ngauss] = info.w[ng:-1:1]
        else
            x, w = gausslegendre(ngauss)
            info.x[1:ngauss] = x
            info.w[1:ngauss] = w
        end

        info.s[1:ngauss] = 1 ./ (sin ∘ acos).(info.x[1:ngauss])

        calc_r!(scatterer, ngauss, info.r, info.dr)
        info.kr1[1:ngauss] = 1.0 ./ (k * info.r[1:ngauss])
        info.kr_s1[1:ngauss] = 1.0 ./ (scatterer.m * k * info.r[1:ngauss])
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.
    rmax = maximum(info.r[1:ngauss])
    kr = k * info.r[1:ngauss]
    kr_s = scatterer.m * kr
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(1.2 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb + 4.0 * ∛tb + 1.2 * √tb - nmax + 5))

    for i in 1:ngauss
        sphericalbesselj!(kr[i], nmax, nnmax1, view(info.jkr, i, :), view(info.djkr, i, :), info.j_tmp)
        sphericalbessely!(kr[i], nmax, view(info.ykr, i, :), view(info.dykr, i, :))
        sphericalbesselj!(kr_s[i], nmax, nnmax2, view(info.jkr_s, i, :), view(info.djkr_s, i, :), info.j_s_tmp)
    end

    info.hkr[1:ngauss, 1:nmax] = complex.(view(info.jkr, 1:ngauss, 1:nmax), view(info.ykr, 1:ngauss, 1:nmax))
    info.dhkr[1:ngauss, 1:nmax] = complex.(view(info.djkr, 1:ngauss, 1:nmax), view(info.dykr, 1:ngauss, 1:nmax))

    info.ngauss = ngauss
    info.nmax = nmax

    return
end

function constant(scatterer::AbstractScatterer, ngauss::Int64, nmax::Int64)
    an = [Float64(n * (n + 1)) for n in 1:nmax]
    ann = [0.5 * √((2n1 + 1) * (2n2 + 1) / (n1 * (n1 + 1) * n2 * (n2 + 1))) for n1 in 1:nmax, n2 in 1:nmax]

    if typeof(scatterer) <: Cylinder
        ng = ngauss ÷ 2
        ng1 = ng ÷ 2
        ng2 = ng - ng1
        x1, w1 = gausslegendre(ng1)
        x2, w2 = gausslegendre(ng2)
        x = zeros(ngauss)
        w = zeros(ngauss)
        xx = -cos(atan(scatterer.d_to_h))
        x[1:ng1] = 0.5(xx + 1.0) * x1 .+ 0.5(xx - 1.0)
        w[1:ng1] = 0.5(xx + 1.0) * w1
        x[(ng1 + 1):ng] = -0.5xx * x2 .+ 0.5xx
        w[(ng1 + 1):ng] = -0.5xx * w2
        x[(ng + 1):ngauss] = -x[ng:-1:1]
        w[(ng + 1):ngauss] = w[ng:-1:1]
    else
        x, w = gausslegendre(ngauss)
    end

    s = 1 ./ (sin ∘ acos).(x)
    ss = s .^ 2

    return x, w, an, ann, s, ss
end

function vary(scatterer::AbstractScatterer, ngauss::Int64, nmax::Int64)
    T = typeof(scatterer.rev)

    r, dr = calc_r(scatterer, ngauss)
    λ = scatterer.λ
    k = 2π / λ
    kr = k * r
    kr1 = 1.0 ./ kr
    kr_s = scatterer.m * kr
    kr_s1 = 1.0 ./ kr_s
    rmax = maximum(r)
    krmax = k * rmax
    tb = max(nmax, krmax * norm(scatterer.m))
    nnmax1 = Int64(floor(1.2 * √(max(krmax, nmax)) + 3.0))
    nnmax2 = Int64(floor(tb + 4.0 * ∛tb + 1.2 * √tb - nmax + 5))

    jkr = zeros(T, ngauss, nmax)
    djkr = zeros(T, ngauss, nmax)
    ykr = zeros(T, ngauss, nmax)
    dykr = zeros(T, ngauss, nmax)
    jkr_s = zeros(Complex{T}, ngauss, nmax)
    djkr_s = zeros(Complex{T}, ngauss, nmax)

    for i in 1:ngauss
        jkr[i, :], djkr[i, :] = sphericalbesselj(kr[i], nmax, nnmax1)
        ykr[i, :], dykr[i, :] = sphericalbessely(kr[i], nmax)
        jkr_s[i, :], djkr_s[i, :] = sphericalbesselj(kr_s[i], nmax, nnmax2)
    end

    return r, dr, kr1, kr_s1, jkr, djkr, ykr, dykr, jkr_s, djkr_s
end

function tmatr0!(scatterer::AbstractScatterer, ngauss::Int64, nmax::Int64)
    T = typeof(scatterer.rev)
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

    for n2 in 1:nmax
        for n1 in 1:nmax
            for i in 1:ngss
                if !(sym && (n1 + n2) % 2 == 1)
                    J12[n1, n2] +=
                        wr2[i] *
                        jkr_s[i, n2] *
                        (
                            dhkr[i, n1] * τ[i, n1] * τ[i, n2] +
                            dr[i] / r[i] * an[n1] * hkr[i, n1] * kr1[i] * d[i, n1] * τ[i, n2]
                        )

                    J21[n1, n2] +=
                        wr2[i] *
                        hkr[i, n1] *
                        (
                            djkr_s[i, n2] * τ[i, n1] * τ[i, n2] +
                            dr[i] / r[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d[i, n2] * τ[i, n1]
                        )

                    RgJ12[n1, n2] +=
                        wr2[i] *
                        jkr_s[i, n2] *
                        (
                            djkr[i, n1] * τ[i, n1] * τ[i, n2] +
                            dr[i] / r[i] * an[n1] * jkr[i, n1] * kr1[i] * d[i, n1] * τ[i, n2]
                        )

                    RgJ21[n1, n2] +=
                        wr2[i] *
                        jkr[i, n1] *
                        (
                            djkr_s[i, n2] * τ[i, n1] * τ[i, n2] +
                            dr[i] / r[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d[i, n2] * τ[i, n1]
                        )
                end
            end
        end
    end

    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann

    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann

    k = 2π / scatterer.λ
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

    Q11 .= kk_s * J21 + kk * J12
    Q22 .= kk_s * J12 + kk * J21

    RgQ11 .= kk_s * RgJ21 + kk * RgJ12
    RgQ22 .= kk_s * RgJ12 + kk * RgJ21

    T = -RgQ * inv(Q)

    return T, Q, RgQ
end

function tmatr!(scatterer::AbstractScatterer, m::Int64, ngauss::Int64, nmax::Int64)
    T = typeof(scatterer.rev)
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

    for n2 in mm:nmax
        for n1 in mm:nmax
            for i in 1:ngss
                if !(sym && (n1 + n2) % 2 == 0)
                    OffsetJ11[n1, n2] +=
                        wr2[i] * hkr[i, n1] * jkr_s[i, n2] * (p[i, n1] * τ[i, n2] + τ[i, n1] * p[i, n2])

                    OffsetJ22[n1, n2] +=
                        wr2[i] * (
                            dhkr[i, n1] * djkr_s[i, n2] * (p[i, n1] * τ[i, n2] + τ[i, n1] * p[i, n2]) +
                            dr[i] / r[i] *
                            (
                                an[n1] * hkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * dhkr[i, n1]
                            ) *
                            p[i, n1] *
                            d[i, n2]
                        )

                    OffsetRgJ11[n1, n2] +=
                        wr2[i] * jkr[i, n1] * jkr_s[i, n2] * (p[i, n1] * τ[i, n2] + τ[i, n1] * p[i, n2])

                    OffsetRgJ22[n1, n2] +=
                        wr2[i] * (
                            djkr[i, n1] * djkr_s[i, n2] * (p[i, n1] * τ[i, n2] + τ[i, n1] * p[i, n2]) +
                            dr[i] / r[i] *
                            (
                                an[n1] * jkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                                an[n2] * jkr_s[i, n2] * kr_s1[i] * djkr[i, n1]
                            ) *
                            p[i, n1] *
                            d[i, n2]
                        )
                end

                if !(sym && (n1 + n2) % 2 == 1)
                    OffsetJ12[n1, n2] +=
                        wr2[i] *
                        jkr_s[i, n2] *
                        (
                            dhkr[i, n1] * (p[i, n1] * p[i, n2] + τ[i, n1] * τ[i, n2]) +
                            dr[i] / r[i] * an[n1] * hkr[i, n1] * kr1[i] * d[i, n1] * τ[i, n2]
                        )

                    OffsetJ21[n1, n2] +=
                        wr2[i] *
                        hkr[i, n1] *
                        (
                            djkr_s[i, n2] * (p[i, n1] * p[i, n2] + τ[i, n1] * τ[i, n2]) +
                            dr[i] / r[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d[i, n2] * τ[i, n1]
                        )

                    OffsetRgJ12[n1, n2] +=
                        wr2[i] *
                        jkr_s[i, n2] *
                        (
                            djkr[i, n1] * (p[i, n1] * p[i, n2] + τ[i, n1] * τ[i, n2]) +
                            dr[i] / r[i] * an[n1] * jkr[i, n1] * kr1[i] * d[i, n1] * τ[i, n2]
                        )

                    OffsetRgJ21[n1, n2] +=
                        wr2[i] *
                        jkr[i, n1] *
                        (
                            djkr_s[i, n2] * (p[i, n1] * p[i, n2] + τ[i, n1] * τ[i, n2]) +
                            dr[i] / r[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d[i, n2] * τ[i, n1]
                        )
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

    k = 2π / scatterer.λ
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

    Q11 .= kk_s * J21 + kk * J12
    Q12 .= kk_s * J11 + kk * J22
    Q21 .= kk_s * J22 + kk * J11
    Q22 .= kk_s * J12 + kk * J21

    RgQ11 .= kk_s * RgJ21 + kk * RgJ12
    RgQ12 .= kk_s * RgJ11 + kk * RgJ22
    RgQ21 .= kk_s * RgJ22 + kk * RgJ11
    RgQ22 .= kk_s * RgJ12 + kk * RgJ21

    T = -RgQ * inv(Q)

    return T, Q, RgQ
end
