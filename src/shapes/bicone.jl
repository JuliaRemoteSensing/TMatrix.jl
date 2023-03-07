@doc raw"""
A bicone scatterer.

Attributes:

- `m`: The complex refractive index.
- `r`: The radius of the cone.
- `h`: The height of one cone.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Bicone{T <: Real, CT <: Number, RV, RM, CV, CM} <: AbstractScatterer{T, CT}
    m::CT
    r::T
    h::T
    λ::T
    info::ScattererInfo{RV, RM, CV, CM}
end

@doc raw"""
Construct a bicone.
"""
function Bicone(T::Type{<:Real} = Float64; r::Real, h::Real, m::Number, λ::Real)
    m = Complex{T}(m)
    r = T(r)
    h = T(h)
    λ = T(λ)
    return Bicone(m, r, h, λ, ScattererInfo(T))
end

has_symmetric_plane(bicone::Bicone) = true

volume_equivalent_radius(bicone::Bicone) = ∛(bicone.r^2 * bicone.h / 2)

function calc_r!(scatterer::Bicone{T},
                 ngauss::Int64,
                 x::AbstractArray,
                 w::AbstractArray,
                 r::AbstractArray,
                 dr::AbstractArray) where {T <: Real}
    theta_split!(scatterer, ngauss, x, w)
    α = atan(scatterer.r / scatterer.h)
    h = scatterer.h
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

function theta_split!(scatterer::Bicone{T}, ngauss::Int64, x::AbstractArray,
                      w::AbstractArray) where {T <: Real}
    ng = ngauss ÷ 2
    x1, w1 = gausslegendre(T, ng)
    @. x[1:ng] = 0.5(x1 - 1)
    @. w[1:ng] = 0.5w1
    @. x[(ng + 1):ngauss] = (-1) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]
end
