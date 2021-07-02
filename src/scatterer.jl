@enum Shape SHAPE_SPHEROID SHAPE_CYLINDER SHAPE_CHEBYSHEV
@enum RadiusType RADIUS_EQUAL_VOLUME RADIUS_EQUAL_AREA RADIUS_MAXIMUM

abstract type AbstractScatterer end

@doc raw"""
A spheroid scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `a_to_c`: The ratio $a/c$ of the horizontal to rotational axes.
"""
struct Spheroid <: AbstractScatterer
    rev::Float64
    m::ComplexF64
    a_to_c::Float64
end

@doc raw"""
A cylinder scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `d_to_h`: The diameter-to-height ratio $D/H$.
"""
struct Cylinder <: AbstractScatterer
    rev::Float64
    m::Complex{Float64}
    d_to_h::Float64
end

@doc raw"""
A Chebyshev scatterer defined by

$$
r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))
$$

in which $T_n(\cos\theta)=\cos n\theta$.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `ε`: The deformation parameter.
- `n`: The degree of the Chebyshev polynomial.
"""
struct Chebyshev <: AbstractScatterer
    rev::Float64
    m::Complex{Float64}
    ε::Float64
    n::Int64
end

@doc raw"""
Scatterer constructor with named parameters.

Parameters:

- `r`: The equivalent radius.
- `radius_type`: The type of the equivalent radius. Possible values are `RADIUS_EQUAL_VOLUME` (default), `RADIUS_EQUAL_AREA` and `RADIUS_MAXIMUM`. All radius types will be transformed into the equivalent volume radius.
- `shape`: The particle shape. Possible values are `SHAPE_SPHEROID` (default), `SHAPE_CYLINDER` and `SHAPE_CHEBYSHEV`.
- `m`: The complex refractive index.
- `axis_ratio`: For spheroids, it is the ratio $a/b$ of the horizontal to rotational axes. For cylinders, it is the diameter-to-length ratio $D/L$.
"""
function Scatterer(; r::Float64, axis_ratio::Float64, shape::Shape=SHAPE_SPHEROID, radius_type::RadiusType=RADIUS_EQUAL_VOLUME, n::Int64=2, refractive_index::ComplexF64=1.0)
    if radius_type == RADIUS_EQUAL_VOLUME
        rev = r
    elseif radius_type == RADIUS_EQUAL_AREA
        # TODO: RADIUS_EQUAL_AREA => RADIUS_EQUAL_VOLUME
        rev = r
    else
        # TODO: RADIUS_MAXIMUM => RADIUS_EQUAL_VOLUME
        rev = r
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
function calc_tmatrix(scatterer::AbstractScatterer, accuracy::Float64=0.001)
end

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
function calc_amplitude(scatterer::AbstractScatterer, tmatrix::Union{Array{Float64,2},Nothing}, ϑ_i::Float64, φ_i::Float64, ϑ_s::Float64, φ_s::Float64)
    if tmatrix === nothing
        tmatrix = calc_tmatrix(scatterer)
    end
end

@doc raw"""
```
calc_r(scatterer::Scatterer, ngauss::Int64)
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer.
"""
function calc_r(scatterer::AbstractScatterer, ngauss::Int64)
    x, _ = gausslegendre(ngauss)

    rev = scatterer.rev
    r = zeros(ngauss)
    dr = zeros(ngauss)

    if typeof(scatterer) == Cylinder
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for cylinders")
        end

        e = scatterer.d_to_h
        h = rev * ∛(2 / (3e^2))
        d = h * e

        @simd for i in 1:ngauss ÷ 2
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

    elseif typeof(scatterer) == Spheroid
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for spheroids")
        end

        # TODO: Need to investigate this equation
        e = scatterer.a_to_c
        a = rev * ∛e

        @simd for i in 1:ngauss ÷ 2
            cosθ = x[i]
            sinθ = √(1.0 - cosθ^2)
            r[i] = a * √(1.0 / (e^2 * cosθ^2 + sinθ^2))
            r[ngauss + 1 - i] = r[i]
            dr[i] = r[i]^3 * cosθ * sinθ * (e^2 - 1.0) / a^2
            dr[ngauss + 1 - i] = -dr[i]
        end
    else
        @assert typeof(scatterer) == Chebyshev
        e = scatterer.ε
        n = scatterer.n
        dn = float(n)

        # TODO: Need to investigate this equation
        a = 1.5e^2 * (4.0dn^2 - 2.0) / (4.0dn^2 - 1.0) + 1.0

        if n % 2 == 0
            a -= 3.0e * (1.0 + 0.25e^2) / (dn^2 - 1.0) + 0.25e^3  / (9.0dn^2 - 1.0)
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
