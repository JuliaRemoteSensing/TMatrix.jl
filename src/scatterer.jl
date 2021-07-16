@enum Shape SHAPE_SPHEROID SHAPE_CYLINDER SHAPE_CHEBYSHEV
@enum RadiusType RADIUS_EQUAL_VOLUME RADIUS_EQUAL_AREA RADIUS_MAXIMUM

abstract type AbstractScatterer end

const CHEBYSHEV_DEFAULT_GAUSSIAN_POINTS = 60

@doc raw"""
A spheroid scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `a_to_c`: The ratio $a/c$ of the horizontal to rotational axes.
"""
struct Spheroid{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    a_to_c::T
end

@doc raw"""
A cylinder scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `d_to_h`: The diameter-to-height ratio $D/H$.
"""
struct Cylinder{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    d_to_h::T
end

@doc raw"""
A Chebyshev scatterer defined by

$r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))$

in which $T_n(\cos\theta)=\cos n\theta$.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `ε`: The deformation parameter.
- `n`: The degree of the Chebyshev polynomial.
"""
struct Chebyshev{T<:Real} <: AbstractScatterer
    rev::T
    m::Complex{T}
    ε::T
    n::Int64
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
) where {T<:Real}
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
            v = 0.0
            @simd for i in 1:ngauss
                θ = acos(x[i])
                nθ = n * θ
                sinθ = sin(θ)
                sinnθ = sin(nθ)
                cosnθ = cos(nθ)
                a = 1.0 + e * cosnθ
                ens = en * sinnθ
                v += w[i] * (sinθ * a + x[i] * ens) * sinθ * a^2
            end
            rv = ∛(0.75v)
            rev = r / (1.0 + e) * rv
        end
    end

    if shape == SHAPE_SPHEROID
        return Spheroid(rev, refractive_index, axis_ratio)
    elseif shape == SHAPE_CYLINDER
        return Cylinder(rev, refractive_index, axis_ratio)
    elseif shape == SHAPE_CHEBYSHEV
        return Chebyshev(rev, refractive_index, axis_ratio, n)
    end
end

@doc raw"""
```
calc_tmatrix(scatterer::Scatterer, accuracy::Float64=0.001)
```

Calculate the T-Matrix of the scatterer.
"""
function calc_tmatrix(scatterer::AbstractScatterer, accuracy::Float64 = 0.001) end

@doc raw"""
```
calc_amplitude(scatterer::Scatterer, tmatrix::Union{Array{Float64,2},Nothing}, ϑ_i::Float64, φ_i::Float64, ϑ_s::Float64, φ_s::Float64)
```

Calculate the amplitude matrix and the phase matrix, given the scatterer and the geometry of the incident and the scattered beam. Use pre-computed T-Matrix when possible.

Parameters:

- `scatterer`: The scatterer.
- `tmatrix`: The pre-computed T-Matrix of the scatterer, or nothing if there is no pre-computation.
- `ϑ_i`: The zenith angle of the incident beam.
- `ϑ_s`: The zenith angle of the scattered beam.
- `φ_i`: The azimuth angle of the indicent beam.
- `φ_s`: The azimuth angle of the scatterer beam.
"""
function calc_amplitude(
    scatterer::AbstractScatterer,
    tmatrix::Union{Array{Float64,2},Nothing},
    ϑ_i::Float64,
    φ_i::Float64,
    ϑ_s::Float64,
    φ_s::Float64,
)
    if tmatrix === nothing
        tmatrix = calc_tmatrix(scatterer)
    end
end

calc_S = calc_amplitude

@doc raw"""

Calculate the phase matrix using the given amplitude matrix $\mathbf{S}$.
"""
function calc_phase(S::Array{Complex{T},2}) where {T<:Real}
    @assert size(S) == (2, 2)

    Z = zeros(T, 4, 4)
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

    return Z
end

calc_Z = calc_phase

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

        # TODO: Need to investigate this equation
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

        # TODO: Need to investigate this equation
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

    return r, dr
end
