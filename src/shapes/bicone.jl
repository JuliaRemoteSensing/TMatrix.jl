@doc raw"""
A bicone scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `r_to_h`: The radius-to-height ratio $R/H$.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Bicone{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    r_to_h::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

function Bicone(T::Type{<:Real} = Float64; r::Real = 1.0, h::Real = 1.0, m::Number = 1.0 + 0.0im, λ::Real = 1.0)
    r = T(r)
    h = T(h)
    m = Complex{T}(m)
    λ = T(λ)
    rev = ∛(r^2 * h / 2)
    return Bicone(rev, m, r / h, λ, ScattererInfo(T))
end

has_symmetric_plane(bicone::Bicone) = true

function calc_r!(
    scatterer::Bicone{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    e = scatterer.r_to_h
    h = rev * ∛(2 / e^2)
    α = atan(e)
    sinα = sin(α)

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = abs(x[i])
        θ = acos(cosθ)
        β = π - α - θ
        sinβ = sin(β)
        cosβ = cos(β)
        r[i] = h / sinβ * sinα
        r[ngauss + 1 - i] = r[i]
        dr[i] = -r[i] * cosβ / sinβ
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function theta_split!(scatterer::Bicone{T}, ngauss::Int64, x::AbstractArray, w::AbstractArray) where {T<:Real}
    ng = ngauss ÷ 2
    x1, w1 = gausslegendre(T, ng)
    @. x[1:ng] = 0.5(x1 - 1)
    @. w[1:ng] = 0.5w1
    @. x[(ng + 1):ngauss] = (-1) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]
end

function calc_r!(
    scatterer::Bicone{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    e = scatterer.r_to_h
    h = rev * ∛(2 / e^2)
    α = atan(e)
    sinα = sin(α)

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = abs(x[i])
        θ = acos(cosθ)
        β = π - α - θ
        sinβ = sin(β)
        cosβ = cos(β)
        r[i] = h / sinβ * sinα
        r[ngauss + 1 - i] = r[i]
        dr[i] = -r[i] * cosβ / sinβ
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function theta_split!(scatterer::Bicone{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike, w::Arblib.ArbVectorLike)
    ng = ngauss ÷ 2
    x1, w1 = gausslegendre(Arb, ng)
    @. x[1:ng] = 0.5(x1 - 1)
    @. w[1:ng] = 0.5w1
    @. x[(ng + 1):ngauss] = (-1) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]
end
