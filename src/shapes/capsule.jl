@doc raw"""
A capsule scatterer, which is a cylinder with half-spherical caps on both ends.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `d`: The diameter of the cylinder base and the hemisphere.
- `h`: The height of the cylinder.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Capsule{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    d::T
    h::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

function Capsule(T::Type{<:Real} = Float64; d::Real = 1.0, h::Real = 1.0, m::Number = 1.0 + 0.0im, λ::Real = 1.0)
    d = T(d)
    h = T(h)
    m = Complex{T}(m)
    λ = T(λ)
    rev = ∛(d^3 / 8 + 3d^2 * h / 16)
    return Capsule(rev, m, d, h, λ, ScattererInfo(T))
end

has_symmetric_plane(capsule::Capsule) = true

@doc raw"""
```
calc_r!(scatterer::Capsule{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given capsule, in place.
"""
function calc_r!(
    scatterer::Capsule{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    d = scatterer.d / 2
    h = scatterer.h / 2

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = abs(x[i])
        sinθ = √(1 - cosθ^2)
        θ = acos(cosθ)
        if h / cosθ < d / sinθ
            sinα = h * sinθ / d
            α = asin(sinα)
            r[i] = d * sin(α + θ) / sinθ
            dr[i] = -h * sinθ - h^2 * sinθ * cosθ / √(d^2 - h^2 * sinθ^2)
        else
            r[i] = d / sinθ
            dr[i] = -d * cosθ / sinθ^2
        end
        r[ngauss + 1 - i] = r[i]
        dr[ngauss + 1 - i] = dr[i]
        dr[i] = -dr[i]
    end
end

function theta_split!(scatterer::Capsule{T}, ngauss::Int64, x::AbstractArray, w::AbstractArray) where {T<:Real}
    ng = ngauss ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1
    x1, w1 = gausslegendre(T, ng1)
    x2, w2 = gausslegendre(T, ng2)
    xx = -cos(atan(scatterer.d / scatterer.h))
    @. x[1:ng1] = 0.5(xx + 1.0) * x1 + 0.5(xx - 1.0)
    @. w[1:ng1] = 0.5(xx + 1.0) * w1
    @. x[(ng1 + 1):ng] = -0.5xx * x2 + 0.5xx
    @. w[(ng1 + 1):ng] = -0.5xx * w2
    @. x[(ng + 1):ngauss] = (-1.0) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]
    return
end
