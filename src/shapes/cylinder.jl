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

function Cylinder(T::Type{<:Real} = Float64; r::Real = 1.0, h::Real = 1.0, m::Number = 1.0 + 0.0im, λ::Real = 1.0)
    r = T(r)
    h = T(h)
    m = Complex{T}(m)
    λ = T(λ)
    rev = ∛(3r^2 * h / 4)
    return Cylinder(rev, m, 2r / h, λ, ScattererInfo(T)) # Axis ratio is defined as D/H for cylinders.
end

has_symmetric_plane(cylinder::Cylinder) = true

@doc raw"""
```
calc_r!(scatterer::AbstractScatterer{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}
```

Calculate $r(\theta)$ and $\frac{\mathrm{d}r}{\mathrm{d}\theta}$ at `ngauss` points for a given scatterer, in place.
"""
function calc_r!(
    scatterer::Cylinder{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
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
end

function theta_split!(scatterer::Cylinder{T}, ngauss::Int64, x::AbstractArray, w::AbstractArray) where {T<:Real}
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
    return
end

function calc_r!(
    scatterer::Cylinder{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    e = scatterer.d_to_h
    h = rev * ∛(2 / (3e^2))
    d = h * e

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

function theta_split!(scatterer::Cylinder{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike, w::Arblib.ArbVectorLike)
    ng = ngauss ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1
    x1, w1 = gausslegendre(Arb, ng1)
    x2, w2 = gausslegendre(Arb, ng2)
    xx = -cos(atan(scatterer.d_to_h))
    x[1:ng1] .= 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
    w[1:ng1] .= 0.5(xx + 1.0) .* w1
    x[(ng1 + 1):ng] .= -0.5xx .* x2 .+ 0.5xx
    w[(ng1 + 1):ng] .= -0.5xx .* w2
    x[(ng + 1):ngauss] .= (-1.0) .* x[ng:-1:1]
    return w[(ng + 1):ngauss] .= w[ng:-1:1]
end
