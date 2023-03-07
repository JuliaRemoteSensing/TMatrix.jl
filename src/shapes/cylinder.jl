@doc raw"""
A cylinder scatterer.

Attributes:

- `m`: The complex refractive index.
- `r`: The radius of the cylinder base.
- `h`: The height of the cylinder.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Cylinder{T <: Real, CT <: Number, RV, RM, CV, CM} <: AbstractScatterer{T, CT}
    m::CT
    r::T
    h::T
    λ::T
    info::ScattererInfo{RV, RM, CV, CM}
end

@doc raw"""
Construct a cylinder. `λ` and `m` must be provided, and the size parameters are considered according to the following priority:

1. `r` and `h`
2. `r` and `r_to_h`
3. `h` and `r_to_h`
4. `rev` and `r_to_h`
5. `rea` and `r_to_h`
6. `rmax` and `r_to_h`

If none of the above is hit, an `ArgumentError` will be thrown.

"""
function Cylinder(T::Type{<:Real} = Float64;
                  r::Real = 0,
                  h::Real = 0,
                  r_to_h::Real = 0,
                  rev::Real = 0,
                  rea::Real = 0,
                  rmax::Real = 0,
                  m::Number,
                  λ::Real)
    m = Complex{T}(m)
    r = T(r)
    h = T(h)
    rev = T(rev)
    rea = T(rea)
    r_to_h = T(r_to_h)
    λ = T(λ)

    if !iszero(r) && !iszero(h)
        nothing
    elseif !iszero(r) && !iszero(r_to_h)
        h = r / r_to_h
    elseif !iszero(h) && !iszero(r_to_h)
        r = h * r_to_h
    elseif !iszero(rev) && !iszero(r_to_h)
        h = rev * ∛(4 / (3r_to_h^2))
        r = h * r_to_h
    elseif !iszero(rea) && !iszero(r_to_h)
        e = 2r_to_h
        ratio = ∛(R_3_2 / e) / √((e + 2) / 2e)
        rev = ratio * rea
        h = rev * ∛(4 / (3r_to_h^2))
        r = h * r_to_h
    elseif !iszero(rmax) && !iszero(r_to_h)
        e = 2r_to_h
        rev = e > 1.0 ? rmax * ∛(R_3_2 / e) : rmax * ∛(R_3_2 * e^2)
        h = rev * ∛(4 / (3r_to_h^2))
        r = h * r_to_h
    else
        throw(ArgumentError("Cannot construct a valid cylinder with the given parameters."))
    end

    return Cylinder(m, r, h, λ, ScattererInfo(T)) # Axis ratio is defined as D/H for cylinders.
end

has_symmetric_plane(cylinder::Cylinder) = true

volume_equivalent_radius(cylinder::Cylinder) = ∛(3cylinder.r^2 * cylinder.h / 4)

@doc raw"""
```
calc_r!(scatterer::AbstractScatterer{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer, in place.
"""
function calc_r!(scatterer::Cylinder{T},
                 ngauss::Int64,
                 x::AbstractArray,
                 w::AbstractArray,
                 r::AbstractArray,
                 dr::AbstractArray) where {T <: Real}
    theta_split!(scatterer, ngauss, x, w)
    h = scatterer.h / 2
    d = scatterer.r

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
end

function theta_split!(scatterer::Cylinder{T}, ngauss::Int64, x::AbstractArray,
                      w::AbstractArray) where {T <: Real}
    ng = ngauss ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1
    x1, w1 = gausslegendre(T, ng1)
    x2, w2 = gausslegendre(T, ng2)
    xx = -cos(atan(2scatterer.r / scatterer.h))
    x[1:ng1] .= 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
    w[1:ng1] .= 0.5(xx + 1.0) .* w1
    x[(ng1 + 1):ng] .= -0.5xx .* x2 .+ 0.5xx
    w[(ng1 + 1):ng] .= -0.5xx .* w2
    x[(ng + 1):ngauss] .= (-1.0) .* x[ng:-1:1]
    w[(ng + 1):ngauss] .= w[ng:-1:1]
    return
end

function calc_r!(scatterer::Cylinder{Arb},
                 ngauss::Int64,
                 x::Arblib.ArbVectorLike,
                 w::Arblib.ArbVectorLike,
                 r::Arblib.ArbVectorLike,
                 dr::Arblib.ArbVectorLike)
    theta_split!(scatterer, ngauss, x, w)
    h = scatterer.h / 2
    d = scatterer.r

    for i in 1:(ngauss ÷ 2)
        cosθ = abs(x[i])
        cos²θ = cosθ^2
        sin²θ = 1 - cos²θ
        sinθ = √sin²θ
        r₁ = h / cosθ
        r₂ = d / sinθ
        if r₁ < r₂
            r[i] = r₁
            dr[i] = h * sinθ / cos²θ
        else
            r[i] = r₂
            dr[i] = -d * cosθ / sin²θ
        end
        r[ngauss + 1 - i] = r[i]
        dr[ngauss + 1 - i] = dr[i]
        dr[i] = -dr[i]
    end
end

function theta_split!(scatterer::Cylinder{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike,
                      w::Arblib.ArbVectorLike)
    ng = ngauss ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1
    x1, w1 = gausslegendre(Arb, ng1)
    x2, w2 = gausslegendre(Arb, ng2)
    xx = -cos(atan(2scatterer.r / scatterer.h))
    x[1:ng1] .= 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
    w[1:ng1] .= 0.5(xx + 1.0) .* w1
    x[(ng1 + 1):ng] .= -0.5xx .* x2 .+ 0.5xx
    w[(ng1 + 1):ng] .= -0.5xx .* w2
    x[(ng + 1):ngauss] .= (-1.0) .* x[ng:-1:1]
    return w[(ng + 1):ngauss] .= w[ng:-1:1]
end
